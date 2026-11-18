#!/usr/bin/env python3
"""
transcript_check.py -- is a given cut site inside EVERY coding isoform, or only
some of them? And in particular: is exon 1 safe to target?

Why this is a separate question from `next_atg`
-----------------------------------------------
The main pipeline already answers "could translation re-initiate at a downstream
in-frame ATG inside this transcript". That is a *within-transcript* escape route.

This script answers the other one: *between* transcripts. If the gene has a
second coding isoform that starts at a different promoter / first exon, then a
frameshift in exon 1 of the reference transcript may simply not be present in
that isoform's mRNA at all, and it keeps making protein. No amount of inDelphi
frameshift score protects you from that.

Both have to be clear before exon 1 is a reasonable target.

Deterministic. It reports isoform geometry and per-guide coverage; which guide
to actually pick is a judgement call -- see the `crispr-primers` skill and
DESIGN_NOTES.

Usage
-----
  python transcript_check.py --gene geneA --chrom chr20 \
         --from out_geneA/geneA_danRer11.csv:3,6,34,15,33,5

  # or just describe the gene, with no guides:
  python transcript_check.py --gene yap1 --chrom chr18 --build danRer11
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

try:
    import certifi
    import ssl
    _CTX = ssl.create_default_context(cafile=certifi.where())
    ssl._create_default_https_context = lambda *a, **k: _CTX  # noqa: SLF001
except Exception:
    pass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chopchop_pipeline import UCSC_API, get_json  # noqa: E402


def coding_len(r: dict) -> int:
    xs = [int(x) for x in r["exonStarts"].rstrip(",").split(",")]
    xe = [int(x) for x in r["exonEnds"].rstrip(",").split(",")]
    return sum(max(0, min(b, r["cdsEnd"]) - max(a, r["cdsStart"]))
               for a, b in zip(xs, xe))


def exons(r: dict) -> list[tuple[int, int]]:
    xs = [int(x) for x in r["exonStarts"].rstrip(",").split(",")]
    xe = [int(x) for x in r["exonEnds"].rstrip(",").split(",")]
    return list(zip(xs, xe))


def fetch_transcripts(build: str, chrom: str, gene: str,
                      lo: int, hi: int) -> list[dict]:
    d = get_json(f"{UCSC_API}/getData/track?genome={build};track=ncbiRefSeq"
                 f";chrom={chrom};start={max(0, lo)};end={hi}")
    rows = d.get("ncbiRefSeq") or []
    if isinstance(rows, dict):
        rows = [r for v in rows.values() for r in v]
    out = [r for r in rows
           if str(r.get("name2", "")).lower() == gene.lower()
           and r.get("cdsStart", 0) < r.get("cdsEnd", 0)]
    return out


def in_cds(r: dict, pos: int) -> bool:
    """Is this genomic position inside the CDS of this transcript (not an intron)?"""
    if not (r["cdsStart"] <= pos < r["cdsEnd"]):
        return False
    return any(a <= pos < b for a, b in exons(r))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Coding-isoform geometry around a locus, and which guides "
                    "hit all of them.")
    ap.add_argument("--gene", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--build", default="danRer11")
    ap.add_argument("--from", dest="source", default=None,
                    metavar="CSV:RANKS",
                    help="pipeline CSV and ranks, e.g. out/yap1_danRer11.csv:5,16")
    ap.add_argument("--span", type=int, default=400000,
                    help="bp around the first transcript found to search")
    ap.add_argument("-o", "--out", default=None, help="write a CSV here too")
    args = ap.parse_args()

    # Locating the gene needs an anchor. The guide CSV already carries exact cut
    # coordinates, so use those; otherwise ask UCSC where the gene is. Sweeping
    # from coordinate 0 would miss any gene further along the chromosome.
    anchor_lo = anchor_hi = None
    if args.source:
        path0 = args.source.rsplit(":", 1)[0]
        cuts = [int(r["cut"]) for r in csv.DictReader(open(path0, newline=""))
                if r.get("cut")]
        if cuts:
            anchor_lo, anchor_hi = min(cuts), max(cuts)
    if anchor_lo is None:
        try:
            hit = get_json(f"{UCSC_API}/search?genome={args.build}"
                           f";search={args.gene}")
            for m in hit.get("positionMatches", []):
                for p in m.get("matches", []):
                    pos = p.get("position", "")
                    if pos.startswith(args.chrom + ":"):
                        a, b = pos.split(":")[1].replace(",", "").split("-")
                        anchor_lo, anchor_hi = int(a), int(b)
                        break
                if anchor_lo is not None:
                    break
        except Exception as e:
            log_err = f"{type(e).__name__}: {e}"
            sys.exit(f"could not locate {args.gene} via UCSC search ({log_err}). "
                     f"Pass --from with a guide CSV instead.")
    if anchor_lo is None:
        sys.exit(f"could not locate {args.gene} on {args.chrom}. "
                 f"Pass --from with a guide CSV so the coordinates are known.")

    seed = fetch_transcripts(args.build, args.chrom, args.gene,
                             anchor_lo - args.span, anchor_hi + args.span)
    if not seed:
        sys.exit(f"no coding ncbiRefSeq transcripts named {args.gene} on "
                 f"{args.chrom} near {anchor_lo}-{anchor_hi} -- check --chrom")
    lo = min(r["txStart"] for r in seed) - 1000
    hi = max(r["txEnd"] for r in seed) + 1000
    txs = fetch_transcripts(args.build, args.chrom, args.gene, lo, hi)
    txs.sort(key=lambda r: (-coding_len(r), 0 if r["name"].startswith("NM_") else 1,
                            r["name"]))

    ref = txs[0]
    strand = ref["strand"]
    ref_len = coding_len(ref)

    print(f"\n{args.gene} ({args.build} {args.chrom}, strand {strand}) -- "
          f"{len(txs)} coding transcript(s)\n")
    print(f"{'transcript':<18}{'CDS nt':>8}{'aa':>7}{'exons':>7}"
          f"{'vs ref':>9}  {'ATG at':>12}  {'stop at':>12}  start")
    print("-" * 96)
    for r in txs:
        cl = coding_len(r)
        # On the minus strand the start codon is at the high end of the CDS.
        atg = r["cdsEnd"] if strand == "-" else r["cdsStart"]
        stop = r["cdsStart"] if strand == "-" else r["cdsEnd"]
        same_start = "same" if atg == (ref["cdsEnd"] if strand == "-"
                                       else ref["cdsStart"]) else "DIFFERENT"
        print(f"{r['name']:<18}{cl:>8}{cl // 3:>7}{len(exons(r)):>7}"
              f"{100.0 * cl / ref_len:>8.0f}%  {atg:>12}  {stop:>12}  {same_start}")

    starts = {(r["cdsEnd"] if strand == "-" else r["cdsStart"]) for r in txs}
    print()
    if len(txs) == 1:
        print("Only one coding isoform annotated. Exon 1 is not at risk from an\n"
              "alternative promoter, so the only re-initiation worry is the\n"
              "within-transcript one the main pipeline already scores (next_atg).")
    elif len(starts) == 1:
        print(f"All {len(txs)} coding isoforms share one start codon. They differ\n"
              "only downstream, so a frameshift in the first coding exon lands in\n"
              "every isoform. Exon 1 is not disqualified by isoform structure.")
    else:
        print(f"WARNING: {len(starts)} distinct start codons across {len(txs)}\n"
              "coding isoforms. An isoform starting downstream would NOT carry a\n"
              "frameshift made in the reference transcript's first coding exon.\n"
              "Target a later exon that is shared by all isoforms instead.")

    if not args.source:
        return

    path, ranks = args.source.rsplit(":", 1)
    wanted = [x.strip() for x in ranks.split(",") if x.strip()]
    rows = {r["rank"]: r for r in csv.DictReader(open(path, newline=""))}
    print(f"\n{'guide':<10}{'cut':>12}{'exon':>6}  isoforms hit        verdict")
    print("-" * 78)
    out_rows = []
    for rk in wanted:
        g = rows.get(rk)
        if not g:
            print(f"rank {rk}: not in {path}")
            continue
        cut = int(g["cut"])
        hit = [r["name"] for r in txs if in_cds(r, cut)]
        n, tot = len(hit), len(txs)
        verdict = ("hits every isoform" if n == tot else
                   f"MISSES {tot - n} isoform(s): "
                   + ", ".join(r["name"] for r in txs if r["name"] not in hit))
        print(f"g{rk:<9}{cut:>12}{g.get('exon', ''):>6}  {n}/{tot:<18} {verdict}")
        out_rows.append(dict(gene=args.gene, rank=rk, cut=cut,
                             exon=g.get("exon", ""), isoforms_hit=n,
                             isoforms_total=tot, verdict=verdict))

    if args.out and out_rows:
        with open(args.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
            w.writeheader()
            w.writerows(out_rows)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
