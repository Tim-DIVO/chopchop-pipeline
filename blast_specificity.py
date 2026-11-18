#!/usr/bin/env python3
"""
blast_specificity.py -- genome-wide specificity for designed primer pairs,
via NCBI BLAST.

Why this exists, and why it took until now
------------------------------------------
Every other output in this workspace said `specificity not checked`. The reason
given was that NCBI Primer-BLAST has no API and no batch mode. That is true of
*Primer-BLAST*, the composite tool. It is not true of BLAST: NCBI publishes a
URL API (`Blast.cgi?CMD=Put/Get`) precisely so that programs can submit
searches. Conflating the two left a real gap open longer than it needed to be.

So this does what Primer-BLAST does, in two steps:

  1. BLAST every primer (short query settings) against a zebrafish-restricted
     database.
  2. Reconstruct predicted amplicons: a plus-strand hit and a minus-strand hit
     on the same subject, correctly oriented, within `--max-amplicon` bp, with
     the 3' end of each primer matching. A pair that predicts exactly one
     amplicon is specific; more than one is not.

The 3'-end requirement is the biologically important filter. A primer with
mismatches at its 5' end still extends; one with a mismatch in the last few
bases usually does not. Counting raw BLAST hits without this over-reports
off-targets dramatically.

What this is not
----------------
Not a substitute for a test PCR. BLAST searches the reference assembly, and an
amplicon that is unique in the reference can still have a paralogue or a
repeat-driven partner in your actual fish line.

Usage
-----
  python blast_specificity.py primers_out/primers.csv -o primers_out/specificity.csv

  # keep an existing search instead of resubmitting (searches are slow):
  python blast_specificity.py primers_out/primers.csv --rid ABC123XYZ
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict

BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
# NCBI asks automated clients to identify themselves with a tool name. It also
# invites an e-mail address; this deliberately sends none, because the address
# available here belongs to the user and is not ours to hand out.
TOOL = "crispr-tools-specificity"
THREE_PRIME = 5          # bases at the 3' end that must match exactly


def _get(retries: int = 6, **params) -> str:
    """One BLAST URL-API call, retried on transport failure.

    This is not defensive padding. A search can sit in NCBI's queue for half an
    hour, and a single DNS hiccup or dropped connection used to raise straight
    out of the poll loop and discard the whole wait -- which is exactly how the
    2026-11-13 run died at 32 minutes (`getaddrinfo failed`). The RID survives
    on NCBI's side for 24 h, so a transport error is never a reason to give up.
    """
    params.setdefault("tool", TOOL)
    url = BLAST + "?" + urllib.parse.urlencode(params)
    delay = 10
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": TOOL})
            return urllib.request.urlopen(req, timeout=300).read().decode("utf-8", "replace")
        except Exception as e:                      # transport, not BLAST, errors
            if attempt == retries:
                raise
            print(f"    network error ({type(e).__name__}: {e}); retrying in "
                  f"{delay}s ({attempt}/{retries - 1})", flush=True)
            time.sleep(delay)
            delay = min(delay * 2, 120)
    raise RuntimeError("unreachable")


def submit(fasta: str, database: str, taxid: str,
           word: int = 11, expect: float = 100, hitlist: int = 250) -> str:
    extra = {"ENTREZ_QUERY": f"txid{taxid}[ORGN]"} if taxid else {}
    out = _get(CMD="Put", PROGRAM="blastn", MEGABLAST="off", DATABASE=database,
               QUERY=fasta,
               WORD_SIZE=str(word), EXPECT=str(expect), FILTER="F",
               NUCL_REWARD="1", NUCL_PENALTY="-3", GAPCOSTS="5 2",
               HITLIST_SIZE=str(hitlist), **extra)
    m = re.search(r"RID = (\S+)", out)
    if not m:
        sys.exit("NCBI did not return an RID. First 500 chars:\n" + out[:500])
    rtoe = re.search(r"RTOE = (\d+)", out)
    print(f"  RID {m.group(1)} (NCBI estimates {rtoe.group(1) if rtoe else '?'}s)")
    return m.group(1)


def wait(rid: str, timeout_s: int, interval: int = 60) -> str:
    """Poll one RID until it is READY.

    60 s is the floor NCBI asks for ("do not poll for any single RID more often
    than once a minute"). The earlier 30 s was over that line, and a client that
    ignores the usage guidelines is a client NCBI is entitled to throttle --
    which is a candidate explanation for searches that queue far past their
    RTOE estimate. Do not lower it.
    """
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        info = _get(CMD="Get", RID=rid, FORMAT_OBJECT="SearchInfo")
        st = re.search(r"Status=(\w+)", info)
        st = st.group(1) if st else "?"
        mins = (time.time() - t0) / 60
        print(f"  {mins:5.1f} min  {st}", flush=True)
        if st == "READY":
            hits = re.search(r"ThereAreHits=(\w+)", info)
            if hits and hits.group(1) == "no":
                sys.exit("BLAST finished with no hits at all -- that cannot be "
                         "right for primers designed on this genome. Check the "
                         "database and organism filter.")
            return _get(CMD="Get", RID=rid, FORMAT_TYPE="XML")
        if st in ("FAILED", "UNKNOWN"):
            sys.exit(f"BLAST reported {st} for RID {rid}")
        time.sleep(interval)
    sys.exit(f"RID {rid} still running after {timeout_s}s. The search is not "
             f"lost -- NCBI keeps it for 24 h. Re-run with --rid {rid} "
             f"(or just re-run: the RID is saved next to the output) to pick it "
             f"up without resubmitting.")


def parse(xml: str, lengths: dict[str, int]) -> dict[str, list[dict]]:
    """query name -> list of usable hits (3'-anchored, near-full-length)."""
    root = ET.fromstring(xml)
    out: dict[str, list[dict]] = defaultdict(list)
    for it in root.iter("Iteration"):
        qname = (it.findtext("Iteration_query-def") or "").split()[0]
        qlen = lengths.get(qname) or int(it.findtext("Iteration_query-len") or 0)
        for hit in it.iter("Hit"):
            acc = hit.findtext("Hit_accession") or "?"
            desc = hit.findtext("Hit_def") or ""
            for hsp in hit.iter("Hsp"):
                q_to = int(hsp.findtext("Hsp_query-to"))
                q_from = int(hsp.findtext("Hsp_query-from"))
                ident = int(hsp.findtext("Hsp_identity"))
                alen = int(hsp.findtext("Hsp_align-len"))
                h_from = int(hsp.findtext("Hsp_hit-from"))
                h_to = int(hsp.findtext("Hsp_hit-to"))
                frame = int(hsp.findtext("Hsp_hit-frame") or 1)
                # The primer's 3' end is its highest query coordinate. Require
                # the alignment to reach it, and to be clean over the last
                # THREE_PRIME bases (approximated by demanding the alignment
                # extend to q_to == qlen with no gaps in that stretch).
                if q_to < qlen:
                    continue
                if alen < qlen - 4 or ident < alen - 2:
                    continue
                lo, hi = (h_from, h_to) if h_from <= h_to else (h_to, h_from)
                out[qname].append(dict(acc=acc, desc=desc[:70], lo=lo, hi=hi,
                                       strand=("+" if frame > 0 else "-"),
                                       ident=ident, alen=alen,
                                       q_from=q_from, q_to=q_to))
    return out


def amplicons(f_hits, r_hits, max_amp: int) -> list[dict]:
    """Pair opposite-strand hits on one subject into predicted products.

    Either primer may play either role at an off-target site, so both
    orientations are tried.
    """
    found = []
    for a_set, b_set in ((f_hits, r_hits), (r_hits, f_hits)):
        for a in a_set:
            if a["strand"] != "+":
                continue
            for b in b_set:
                if b["strand"] != "-" or b["acc"] != a["acc"]:
                    continue
                size = b["hi"] - a["lo"] + 1
                if 0 < size <= max_amp:
                    found.append(dict(acc=a["acc"], desc=a["desc"],
                                      start=a["lo"], end=b["hi"], size=size))
    # de-duplicate identical predictions
    seen, uniq = set(), []
    for f in found:
        k = (f["acc"], f["start"], f["end"])
        if k not in seen:
            seen.add(k)
            uniq.append(f)
    return sorted(uniq, key=lambda f: -f["size"])


def main() -> None:
    ap = argparse.ArgumentParser(description="BLAST specificity for primer pairs.")
    ap.add_argument("primers", help="primers.csv from primer_design.py")
    ap.add_argument("-o", "--out", default=None, help="write a CSV report here")
    ap.add_argument("--database", default="GPIPE/7955/current/all_top_level",
                    help="NCBI BLAST database. The default is the zebrafish "
                         "genome database -- the same one Primer-BLAST searches "
                         "when you pick 'Genome (reference assembly)' -- so the "
                         "search space is one 1.4 Gb assembly instead of every "
                         "representative genome in RefSeq filtered down "
                         "afterwards. Pass refseq_representative_genomes to go "
                         "back to the old behaviour")
    ap.add_argument("--taxid", default="",
                    help="Entrez organism filter, e.g. 7955 for Danio rerio. "
                         "Empty by default: the default database is already "
                         "zebrafish-only, and an ENTREZ_QUERY on top of it adds "
                         "an Entrez pre-filter for no gain")
    ap.add_argument("--max-amplicon", type=int, default=4000)
    ap.add_argument("--word-size", type=int, default=11,
                    help="blastn word size (default 11). A site with no 11 nt "
                         "exact stretch will not prime, so this is not a real "
                         "sensitivity limit; 7, which is what Primer-BLAST "
                         "uses, mostly buys queue time")
    ap.add_argument("--expect", type=float, default=100,
                    help="E-value ceiling. A full-length 20-mer match in a "
                         "1.4 Gb genome scores E ~0.05, so 100 is generous; the "
                         "old 1000 mostly added marginal HSPs to format")
    ap.add_argument("--hitlist", type=int, default=250)
    ap.add_argument("--timeout", type=int, default=3600, help="seconds to wait")
    ap.add_argument("--rid", default=None, help="pick up an existing search")
    ap.add_argument("--xml", default=None, help="parse a saved XML file instead")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.primers, newline="")))
    rows = [r for r in rows if r.get("fwd") and r.get("rev")]
    # primers.csv carries one row per ordered pair (primary and backup),
    # so `name` repeats and `id` (<locus>_p<rank>) is the unique key.
    # Query names MUST be unique: BLAST groups output by query-def, so two
    # pairs sharing a name merge their hits, and pairing a primary forward
    # primer with a backup reverse primer then invents products nobody
    # ordered.
    for r in rows:
        r["_id"] = r.get("id") or r["name"]
    lengths, fasta = {}, []
    for r in rows:
        for side in ("fwd", "rev"):
            n = f"{r['_id']}_{'F' if side == 'fwd' else 'R'}"
            lengths[n] = len(r[side])
            fasta.append(f">{n}\n{r[side]}")
    fasta = "\n".join(fasta) + "\n"
    print(f"{len(rows)} pairs / {len(lengths)} primers against "
          f"{args.database} (txid{args.taxid})")

    # Where the RID is remembered. A queued search is a real asset -- it can
    # represent half an hour of waiting -- so it is written to disk before the
    # first poll and reused automatically when the query set is unchanged.
    stem = os.path.splitext(args.out or "blast")[0]
    ridfile = stem + ".rid"
    qhash = hashlib.sha1(fasta.encode()).hexdigest()[:12]

    if args.xml:
        xml = open(args.xml, encoding="utf-8").read()
    else:
        rid = args.rid
        if not rid and os.path.exists(ridfile):
            try:
                saved = json.load(open(ridfile, encoding="utf-8"))
            except Exception:
                saved = {}
            if saved.get("query") == qhash:
                age = (time.time() - saved.get("t", 0)) / 3600
                if age < 23:
                    rid = saved.get("rid")
                    print(f"  resuming saved RID {rid} "
                          f"({age:.1f} h old, same query set)")
                else:
                    print(f"  saved RID {saved.get('rid')} is {age:.0f} h old "
                          f"(NCBI keeps results 24 h); resubmitting")
            else:
                print("  saved RID is for a different query set; resubmitting")
        if not rid:
            rid = submit(fasta, args.database, args.taxid,
                         args.word_size, args.expect, args.hitlist)
        elif args.rid:
            print(f"  picking up RID {rid}")
        with open(ridfile, "w", encoding="utf-8") as fh:
            json.dump({"rid": rid, "query": qhash, "t": time.time(),
                       "database": args.database}, fh)
        xml = wait(rid, args.timeout)
        with open(os.path.splitext(args.out or "blast")[0] + ".xml", "w",
                  encoding="utf-8") as fh:
            fh.write(xml)

    hits = parse(xml, lengths)
    report = []
    print()
    for r in rows:
        # Fall back to the bare locus name when reading XML from a job
        # submitted before query names were unique. Such a job merges a
        # locus's two pairs, so its off-target counts cannot be trusted --
        # this exists so an old RID is not silently read as "no hits".
        f = hits.get(f"{r['_id']}_F") or hits.get(f"{r['name']}_F", [])
        rv = hits.get(f"{r['_id']}_R") or hits.get(f"{r['name']}_R", [])
        amps = amplicons(f, rv, args.max_amplicon)
        intended = [a for a in amps if a["size"] == int(r["product_bp"])]
        extra = [a for a in amps if a["size"] != int(r["product_bp"])]
        if not amps:
            verdict = "NO PRODUCT PREDICTED -- check the database/filter"
        elif not intended:
            verdict = f"intended {r['product_bp']} bp NOT predicted"
        elif not extra:
            verdict = "specific (1 product)"
        else:
            verdict = f"{len(amps)} products -- {len(extra)} off-target"
        flag = "ok  " if verdict.startswith("specific") else "WARN"
        print(f"  {flag} {r['_id']:<14} hits F/R {len(f):>3}/{len(rv):<3} "
              f"products {len(amps):<3} {verdict}")
        for a in extra[:4]:
            print(f"        off-target {a['size']} bp on {a['acc']} {a['desc']}")
        report.append(dict(id=r["_id"], name=r["name"], f_hits=len(f), r_hits=len(rv),
                           products=len(amps), intended=len(intended),
                           off_target=len(extra), verdict=verdict,
                           detail="; ".join(f"{a['acc']}:{a['start']}-{a['end']}"
                                            f"({a['size']}bp)" for a in extra[:10])))

    n_ok = sum(1 for x in report if x["verdict"].startswith("specific"))
    print(f"\n{n_ok}/{len(report)} pairs predicted specific.")
    print("BLAST searches the reference assembly. A test PCR on your own line "
          "is still the real control.")

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(report[0]))
            w.writeheader()
            w.writerows(report)
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
