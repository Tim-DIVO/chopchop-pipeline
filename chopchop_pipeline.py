#!/usr/bin/env python3
"""
chopchop_pipeline.py -- end-to-end CRISPR guide selection for zebrafish F0 knockouts.

Follows Kroll & Rihel's protocols:
  gRNA selection : dx.doi.org/10.17504/protocols.io.81wgb6r5qlpk/v1
  F0 knockout    : dx.doi.org/10.17504/protocols.io.5qpvo52wdl4o/v3
  Method paper   : 10.7554/eLife.59683

What it does
  1. Submits a CHOPCHOP job (or reuses one you pass with --job-url) and waits for it.
  2. Pulls the ranked guide table, inDelphi repair predictions, and primer pairs.
  3. Annotates every guide with exon number, symmetric/asymmetric, and distance to
     the nearest splice boundary, using the RefSeq gene model from the UCSC API.
  4. Lifts each cut site to Zv9/danRer7 with a UCSC chain file and queries SNPfisher
     (one batched request) for strain SNPs across the protospacer+PAM.
  5. Writes an annotated CSV and a self-contained interactive HTML report.

IMPORTANT -- undocumented endpoints
  The per-guide endpoints (`<rank>_repStats.json`, `primer_<rank>.json`) and the job
  submission payload are NOT a supported public API. They were derived from CHOPCHOP's
  own client-side JavaScript. They work as of 2026-11, but may break without notice.
  This script fails loudly rather than emitting blank columns -- that is deliberate.

Example
  python chopchop_pipeline.py --gene yap1 --build danRer11 -o out/
  python chopchop_pipeline.py --job-url https://chopchop.cbu.uib.no/results/<id>/ -o out/
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import webbrowser

CHOPCHOP = "https://chopchop.cbu.uib.no"
UCSC_API = "https://api.genome.ucsc.edu"
UCSC_DL = "https://hgdownload.soe.ucsc.edu/goldenPath"
SNPFISHER = "https://snpfisher.nichd.nih.gov/snpfisher"
UA = "chopchop-pipeline/1.0 (academic guide design; contact: lab user)"

# Builds for which a SNPfisher (Zv9/danRer7) liftover is meaningful.
SNP_CHAIN = {"danRer11": ("danRer11ToDanRer7.over.chain.gz", "danRer7"),
             "danRer10": ("danRer10ToDanRer7.over.chain.gz", "danRer7")}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


SNP_STRAINS = {"fli": "FLI_FQ", "tl": "TL_FQ", "wik": "WIK_FQ"}


def _fq(v) -> float:
    """SNPfisher writes allele frequencies as '0', '1', '0.57' or '\\N' (= no call)."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def mm_int(v: str) -> int:
    """CHOPCHOP writes '>=121' when it hits its off-target cap -- treat as the number."""
    s = str(v).strip().lstrip(">=").strip()
    try:
        return int(s)
    except ValueError:
        raise PipelineError(f"cannot parse off-target count {v!r}")


# ----------------------------------------------------------------------------- errors
class PipelineError(RuntimeError):
    pass


def fail(msg: str, hint: str | None = None) -> "never":
    sys.stderr.write(f"\nFATAL: {msg}\n")
    if hint:
        sys.stderr.write(f"  hint: {hint}\n")
    sys.exit(1)


def _open(req, timeout=120):
    try:
        return urllib.request.urlopen(req, timeout=timeout)
    except urllib.error.HTTPError as e:
        raise PipelineError(f"HTTP {e.code} for {getattr(req,'full_url',req)}") from e
    except urllib.error.URLError as e:
        raise PipelineError(f"network error for {getattr(req,'full_url',req)}: {e.reason}") from e


def get_text(url: str, timeout=120) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    body = _open(req, timeout).read().decode("utf-8", "replace")
    if body.lstrip().lower().startswith("<!doctype html"):
        raise PipelineError(f"expected data but got an HTML page from {url}")
    return body


def get_json(url: str, timeout=120):
    txt = get_text(url, timeout)
    try:
        return json.loads(txt)
    except json.JSONDecodeError as e:
        raise PipelineError(f"invalid JSON from {url}: {e}") from e


def post_json(url: str, payload: dict, timeout=300):
    req = urllib.request.Request(
        url, data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "User-Agent": UA})
    return json.loads(_open(req, timeout).read().decode())


# ------------------------------------------------------------------------ chopchop job
def build_opts(build: str, scoring: str, pam: str, product_min: int, product_max: int) -> list[str]:
    """Mirror the option array CHOPCHOP's own web form sends. -P/-3 produce primers;
    -repairPredictions produces the inDelphi repStats this pipeline depends on."""
    primer3 = (f"PRODUCT_SIZE_MIN={product_min},PRODUCT_SIZE_MAX={product_max},"
               "PRIMER_MIN_SIZE=18,PRIMER_MAX_SIZE=25,PRIMER_OPT_SIZE=22,"
               "PRIMER_MIN_TM=57,PRIMER_MAX_TM=63,PRIMER_OPT_TM=60")
    return ["-J", "-BED", "-GenBank", "-G", build,
            "-filterGCmin", "10", "-filterGCmax", "90",
            "-t", "CODING", "-n", "N", "-R", "4",
            "-3", primer3, "-P", "-A", str(product_max), "-DF", "300",
            "-a", "20", "-T", "1", "-g", "20",
            "-scoringMethod", scoring, "-f", "NN", "-v", "3",
            "-repairPredictions", "mESC", "-M", pam, "-BB", "AGGCTAGTCCGT"]


def submit(gene: str, opts: list[str]) -> str:
    payload = {"opts": opts, "fastaInput": "", "geneInput": gene,
               "isIsoform": False, "forSelect": "knock-out"}
    out = post_json(CHOPCHOP + "/", payload)
    jid = out.get("jobId")
    if not jid:
        raise PipelineError(f"submission returned no jobId: {out}")
    return jid


def wait_for(job: str, timeout_s=900, interval=8) -> str:
    """Poll results.tsv until the job finishes. CHOPCHOP 404s while still running."""
    url = f"{CHOPCHOP}/results/{job}/results.tsv"
    t0 = time.time()
    while time.time() - t0 < timeout_s:
        try:
            return get_text(url, timeout=60)
        except PipelineError as e:
            if "HTTP 404" not in str(e):
                raise
            time.sleep(interval)
    raise PipelineError(f"job {job} not finished after {timeout_s}s")


def parse_results(tsv: str) -> list[dict]:
    rows = list(csv.DictReader(tsv.splitlines(), delimiter="\t"))
    if not rows:
        raise PipelineError("results.tsv is empty -- CHOPCHOP found no target sites")
    need = {"Rank", "Target sequence", "Genomic location", "Strand", "Efficiency",
            "MM0", "MM1", "MM2", "MM3", "GC content (%)", "Self-complementarity"}
    missing = need - set(rows[0])
    if missing:
        raise PipelineError(f"results.tsv missing expected columns: {sorted(missing)}",)
    return rows


# -------------------------------------------------------------------------- gene model
def ucsc_seq(build: str, chrom: str, start0: int, end0: int) -> str:
    d = get_json(f"{UCSC_API}/getData/sequence?genome={build};chrom={chrom}"
                 f";start={start0};end={end0}")
    if "dna" not in d:
        raise PipelineError(f"UCSC returned no sequence for {chrom}:{start0}-{end0}: {d}")
    return d["dna"].upper()


def gene_model(build: str, chrom: str, lo: int, hi: int, gene: str | None):
    """Return (transcript_name, strand, exons[(s,e)] in transcription order, cds0, cds1)
    for the reference transcript: the one with the longest CDS, per Kroll."""
    d = get_json(f"{UCSC_API}/getData/track?genome={build};track=ncbiRefSeq"
                 f";chrom={chrom};start={max(0,lo)};end={hi}")
    rows = d.get("ncbiRefSeq")
    if rows is None:
        raise PipelineError("UCSC ncbiRefSeq track unavailable for this build/region",
                            "check that --build names a UCSC assembly (e.g. danRer11)")
    if isinstance(rows, dict):
        rows = [r for v in rows.values() for r in v]
    cand = [r for r in rows if r.get("cdsStart", 0) < r.get("cdsEnd", 0)]
    if gene:
        named = [r for r in cand if str(r.get("name2", "")).lower() == gene.lower()]
        if named:
            cand = named
    if not cand:
        return None
    def _coding_len(r):
        """True protein-coding length: sum of each exon's overlap with the CDS.
        NOT cdsEnd - cdsStart, which is the genomic span and includes introns."""
        xs = [int(x) for x in r["exonStarts"].rstrip(",").split(",")]
        xe = [int(x) for x in r["exonEnds"].rstrip(",").split(",")]
        return sum(max(0, min(b, r["cdsEnd"]) - max(a, r["cdsStart"]))
                   for a, b in zip(xs, xe))

    # Kroll: "Choose the longest protein-coding transcript as reference" -- longest
    # PROTEIN, so rank by coding length. Ties are common (a curated NM_ and a
    # predicted XM_ describing the same thing), so prefer the curated record, then
    # fall back to the name for determinism.
    def _rank(r):
        curated = 0 if r["name"].startswith(("NM_", "NR_")) else 1
        return (-_coding_len(r), curated, r["name"])

    best = min(cand, key=_rank)
    st = [int(x) for x in best["exonStarts"].rstrip(",").split(",")]
    en = [int(x) for x in best["exonEnds"].rstrip(",").split(",")]
    ex = sorted(zip(st, en), reverse=(best["strand"] == "-"))
    return (best["name"], best["strand"], ex, best["cdsStart"], best["cdsEnd"],
            _coding_len(best))


def annotate_exons(guides: list[dict], model) -> None:
    """Assign each guide's CUT SITE (not its protospacer) to an exon.

    A guide can overlap an exon while cutting a few bp past the splice site, in the
    intron -- an indel there usually will not shift the reading frame. Those are
    reported with region='intron' and a NEGATIVE bp_to_splice (how far past the
    boundary the break falls), so any filter of bp_to_splice >= 0 drops them.
    """
    if model is None:
        for g in guides:
            g.update(exon="", exon_total="", symmetry="", bp_to_splice="",
                     in_cds="", region="")
        return
    _, strand, exons, cds0, cds1, _cdslen = model
    total = len(exons)
    for g in guides:
        cut0 = g["cut"] - 1                      # 0-based
        g.update(exon="", exon_total=total, symmetry="", bp_to_splice="",
                 in_cds="", region="")
        inside = None
        for i, (a, b) in enumerate(exons, 1):    # a,b are 0-based half-open
            if a <= cut0 < b:
                inside = (i, a, b)
                break
        if inside:
            i, a, b = inside
            coding = max(0, min(b, cds1) - max(a, cds0))
            g["exon"] = i
            g["symmetry"] = "asymmetric" if coding % 3 else "symmetric"
            g["bp_to_splice"] = min(cut0 - a, b - 1 - cut0)
            g["in_cds"] = "yes" if cds0 <= cut0 < cds1 else "no"
            g["region"] = "exon"
            continue
        # outside every exon: record the nearest one and how far past it we are
        best_i, best_d = None, None
        for i, (a, b) in enumerate(exons, 1):
            d = 0 if a <= cut0 < b else min(abs(cut0 - a), abs(cut0 - (b - 1)))
            if best_d is None or d < best_d:
                best_i, best_d = i, d
        g["exon"] = best_i
        g["bp_to_splice"] = -best_d
        g["in_cds"] = "no"
        g["region"] = "intron" if best_d <= 2000 else "outside"


# --------------------------------------------------------------------------- liftover
def ensure_chain(build: str, cache_dir: str) -> str | None:
    if build not in SNP_CHAIN:
        return None
    name, _ = SNP_CHAIN[build]
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, name)
    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        url = f"{UCSC_DL}/{build}/liftOver/{name}"
        sys.stderr.write(f"  downloading chain {name} ...\n")
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        data = _open(req, 300).read()
        if len(data) < 1000:
            raise PipelineError(f"chain file from {url} looks truncated ({len(data)} B)")
        with open(path, "wb") as fh:
            fh.write(data)
    return path


def load_chain(path: str, want_t: str):
    blocks, tName, active = [], None, False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("chain"):
                p = line.split()
                tName, tStrand, tPos = p[2], p[4], int(p[5])
                qName, qSize, qStrand, qPos = p[7], int(p[8]), p[9], int(p[10])
                active = (tName == want_t and tStrand == "+")
            elif line and tName:
                p = line.split()
                size = int(p[0])
                if active:
                    blocks.append((tPos, tPos + size, qPos, qName, qStrand, qSize))
                if len(p) == 3:
                    tPos += size + int(p[1]); qPos += size + int(p[2])
                else:
                    tName, active = None, False
    blocks.sort()
    return blocks


def lift(blocks, pos0: int):
    lo, hi = 0, len(blocks) - 1
    while lo <= hi:
        m = (lo + hi) // 2
        ts, te, qs, qn, qstrand, qsize = blocks[m]
        if pos0 < ts:
            hi = m - 1
        elif pos0 >= te:
            lo = m + 1
        else:
            q = qs + (pos0 - ts)
            if qstrand == "-":
                q = qsize - 1 - q
            return qn, q
    return None


# -------------------------------------------------------------------------- snpfisher
def snpfisher_fields(cache_dir: str) -> dict:
    """Scrape the form's own defaults, then set the frequency comparators to '>='.
    The page ships them as '= 0.0', which matches only SNPs at exactly zero frequency
    in all three strains -- i.e. it silently returns nothing. This is the one real trap."""
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, "snpfisher_form.html")
    if not os.path.exists(path):
        html_txt = urllib.request.urlopen(
            urllib.request.Request(f"{SNPFISHER}/snpfisher.html", headers={"User-Agent": UA}),
            timeout=90).read().decode("utf-8", "replace")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(html_txt)
    else:
        html_txt = open(path, encoding="utf-8", errors="replace").read()
    i, j = html_txt.find("<form"), html_txt.find("</form>")
    if i < 0 or j < 0:
        raise PipelineError("could not locate the SNPfisher form -- page layout changed")
    frm = html_txt[i:j]
    fields: dict[str, str] = {}
    for m in re.finditer(r"<input([^>]*)>", frm, re.I):
        a = m.group(1)
        nm = re.search(r"name=[\"']?(\w+)", a)
        if not nm:
            continue
        nm = nm.group(1)
        tm = re.search(r"type=[\"']?(\w+)", a)
        typ = (tm.group(1) if tm else "text").lower()
        vm = re.search(r"value=[\"']([^\"']*)[\"']", a)
        val = vm.group(1) if vm else ""
        if typ == "radio":
            if "checked" in a.lower():
                fields[nm] = val
            else:
                fields.setdefault(nm, "")
        else:
            fields[nm] = val
    for s in ("FLI", "TL", "WIK"):
        fields[f"{s}filter1"] = ">="
        fields[f"{s}filter2"] = "0.0"
        fields[f"{s}filter3"] = "0"
    if "query" not in fields:
        fields["query"] = ""
    return fields


def snpfisher_batch(fields: dict, windows: list[str]) -> dict[str, list[dict]]:
    """One POST for all windows. Output is delimited by '#chrN from X to Y' headers."""
    if not windows:
        return {}
    f = dict(fields)
    f["query"] = "\n".join(windows)
    data = urllib.parse.urlencode(f).encode()
    req = urllib.request.Request(f"{SNPFISHER}/snpTextOracle.jsp", data=data,
                                 headers={"Content-Type": "application/x-www-form-urlencoded",
                                          "User-Agent": UA})
    txt = _open(req, 300).read().decode("utf-8", "replace")
    txt = re.sub(r"<[^>]+>", "", txt)
    if "FLI_FQ" not in txt:
        raise PipelineError("SNPfisher response has no recognisable header",
                            "the form fields or endpoint may have changed")
    cols = ["Chr", "Pos", "Ref", "Mut", "FLI_GT", "FLI_DP", "FLI_AD", "FLI_FQ",
            "TL_GT", "TL_DP", "TL_AD", "TL_FQ", "WIK_GT", "WIK_DP", "WIK_AD",
            "WIK_FQ", "RE_Change"]
    out: dict[str, list[dict]] = {}
    key = None
    for line in txt.splitlines():
        line = line.rstrip()
        if not line.strip():
            continue
        m = re.match(r"#(\w+) from (\d+) to (\d+)", line)
        if m:
            key = f"{m.group(1)}:{m.group(2)}-{m.group(3)}"
            out[key] = []
            continue
        if line.startswith("#") or key is None:
            continue
        p = line.split("\t")
        if len(p) < 16 or p[0] == "Chr":
            continue
        out[key].append(dict(zip(cols, p)))
    return out


# ------------------------------------------------------------------------------- main
def run(args) -> None:
    # Cache lives next to the script, not in the output dir, so the ~2 MB liftOver
    # chain and the SNPfisher form are downloaded once per machine rather than once
    # per gene.
    cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".cache")
    os.makedirs(args.out, exist_ok=True)

    # 1. job -------------------------------------------------------------------
    if args.job_url:
        m = re.search(r"/results/([0-9a-f-]{8,})", args.job_url)
        if not m:
            fail(f"could not parse a job id out of {args.job_url!r}",
                 "expected .../results/<uuid>/")
        job = m.group(1)
        sys.stderr.write(f"[1/6] using existing job {job}\n")
    else:
        if not args.gene:
            fail("need --gene (or --job-url)")
        opts = build_opts(args.build, args.scoring, args.pam,
                          args.product_min, args.product_max)
        sys.stderr.write(f"[1/6] submitting {args.gene} / {args.build} / {args.scoring} ...\n")
        job = submit(args.gene, opts)
        sys.stderr.write(f"      jobId {job}\n")

    tsv = wait_for(job, args.timeout)
    guides = parse_results(tsv)
    sys.stderr.write(f"      {len(guides)} guides\n")

    try:
        query = get_json(f"{CHOPCHOP}/results/{job}/query.json")
    except PipelineError:
        query = {}
    gene = args.gene or query.get("geneInput") or "gene"
    build = args.build
    if query.get("opts") and "-G" in query["opts"]:
        build = query["opts"][query["opts"].index("-G") + 1]

    # 2. normalise coordinates -------------------------------------------------
    sys.stderr.write("[2/6] resolving coordinates against the genome ...\n")
    chroms = {g["Genomic location"].split(":")[0] for g in guides}
    if len(chroms) != 1:
        sys.stderr.write(f"      note: guides span {len(chroms)} sequences: {sorted(chroms)}\n")
    main_chrom = sorted(chroms, key=lambda c: -sum(
        1 for g in guides if g["Genomic location"].startswith(c + ":")))[0]
    pos = [int(g["Genomic location"].split(":")[1]) for g in guides
           if g["Genomic location"].startswith(main_chrom + ":")]
    lo, hi = min(pos) - 500, max(pos) + 500
    ref = ucsc_seq(build, main_chrom, lo - 1, hi)

    kept = []
    for g in guides:
        ch, p = g["Genomic location"].split(":")
        p = int(p)
        if ch != main_chrom or not (lo <= p <= hi - 23):
            continue
        seq = g["Target sequence"].upper()
        win = ref[p - lo: p - lo + len(seq)]
        if win == seq:
            strand = "+"
        elif rc(win) == seq:
            strand = "-"
        else:
            raise PipelineError(
                f"guide {g['Rank']} ({seq}) does not match {build} at {ch}:{p}",
                "is --build the same assembly the job was run on?")
        # Cas9 blunt cut: 3 bp 5' of the PAM. Report the base left of the break.
        cut = p + 16 if strand == "+" else p + 5
        kept.append(dict(rank=int(g["Rank"]), target=seq, proto=seq[:-3], pam=seq[-3:],
                         chrom=ch, start=p, strand=strand, cut=cut,
                         efficiency=float(g["Efficiency"]), gc=int(g["GC content (%)"]),
                         selfcompl=int(g["Self-complementarity"]),
                         mm0=mm_int(g["MM0"]), mm1=mm_int(g["MM1"]),
                         mm2=mm_int(g["MM2"]), mm3=mm_int(g["MM3"]),
                         mm_capped="yes" if ">=" in "".join(
                             str(g[k]) for k in ("MM0", "MM1", "MM2", "MM3")) else ""))
    if not kept:
        raise PipelineError("no guides could be matched to the genome")
    guides = sorted(kept, key=lambda g: g["rank"])
    sys.stderr.write(f"      {len(guides)} guides verified against {build}\n")

    # 3. exon annotation -------------------------------------------------------
    sys.stderr.write("[3/6] annotating exons ...\n")
    model = gene_model(build, main_chrom, lo, hi, gene)
    annotate_exons(guides, model)
    if model:
        sys.stderr.write(f"      reference transcript {model[0]} "
                         f"({len(model[2])} exons, {model[1]} strand, "
                         f"{model[5]} bp CDS = {model[5]//3-1} aa, longest protein)\n")
    else:
        sys.stderr.write("      no coding RefSeq transcript found; exon columns blank\n")

    # 4. inDelphi + primers ----------------------------------------------------
    sys.stderr.write(f"[4/6] fetching inDelphi + primers for {len(guides)} guides ...\n")
    miss_rep = miss_prim = got_primers = 0
    for i, g in enumerate(guides):
        try:
            st = get_json(f"{CHOPCHOP}/results/{job}/{g['rank']}_repStats.json", timeout=60)
            g["frameshift"] = round(float(st["Frameshift frequency"]), 2)
            g["precision"] = round(float(st["Precision"]), 3)
            g["mh_del"] = round(float(st["MH del frequency"]), 2)
            g["exp_indel"] = round(float(st["Expected indel length"]), 2)
        except (PipelineError, KeyError, ValueError):
            g.update(frameshift="", precision="", mh_del="", exp_indel="")
            miss_rep += 1
        try:
            pr = get_json(f"{CHOPCHOP}/results/{job}/primer_{g['rank']}.json", timeout=60)
            if pr and not all(isinstance(r, list) and len(r) >= 15 for r in pr):
                raise PipelineError(
                    f"primer_{g['rank']}.json has an unexpected row shape "
                    f"(got {[len(r) for r in pr[:3]]}, expected >=15 fields)",
                    "CHOPCHOP's primer JSON layout changed; update the field indices")
            # layout: [chrom, fS, fE, rS, rE, idx, strand, fwd, rev, ?, ?, ?, tmF, tmR, size]
            g["primers"] = [dict(fwd=r[7], rev=r[8], tm_f=r[12], tm_r=r[13], size=r[14],
                                 fwd_start=r[1], rev_end=r[4]) for r in pr]
            if pr:
                got_primers += 1
        except PipelineError:
            g["primers"] = []
            miss_prim += 1
        if i % 25 == 0 and i:
            sys.stderr.write(f"      {i}/{len(guides)}\n")
        time.sleep(args.delay)
    if miss_rep == len(guides):
        raise PipelineError("no inDelphi data for ANY guide",
                            "the job may have run without -repairPredictions, or the "
                            "repStats endpoint changed")
    if miss_rep:
        sys.stderr.write(f"      warning: {miss_rep} guides lack inDelphi data\n")
    if got_primers == 0:
        if miss_prim == len(guides):
            sys.stderr.write("      warning: primer endpoint returned nothing for any guide "
                             "(job may have run without -P)\n")
        else:
            raise PipelineError(
                "primer files were fetched but zero primer pairs were parsed",
                "the primer JSON layout has probably changed -- check the field indices "
                "in the primer block rather than shipping a CSV with empty primer columns")
    else:
        sys.stderr.write(f"      primers for {got_primers}/{len(guides)} guides\n")

    # 5. liftover + SNPfisher --------------------------------------------------
    sys.stderr.write("[5/6] liftover to Zv9 + SNPfisher ...\n")
    for g in guides:
        g.update(zv9="", snps="", snp_detail=[],
                 **{f"snp_{k}": "" for k in SNP_STRAINS})
    snp_fields: dict = {}
    chain_path = ensure_chain(build, cache)
    if chain_path is None:
        sys.stderr.write(f"      no Zv9 chain for {build}; skipping SNP check\n")
    else:
        blocks = load_chain(chain_path, main_chrom)
        if not blocks:
            raise PipelineError(f"chain file has no blocks for {main_chrom}")
        _, zbuild = SNP_CHAIN[build]

        # Lift every guide (cheap, local) so the HTML's live SNP button works on any
        # row -- then validate all of them against a single fetched Zv9 span.
        lifted = {}
        for g in guides:
            a = lift(blocks, g["start"] - 1)
            b = lift(blocks, g["start"] + 21)
            if a and b and a[0] == b[0]:
                lifted[g["rank"]] = (a[0], min(a[1], b[1]), max(a[1], b[1]) + 1)
            else:
                g["zv9"] = "liftover-failed"
        spans: dict[str, list[int]] = {}
        for qn, s0, e0 in lifted.values():
            s = spans.setdefault(qn, [s0, e0])
            s[0] = min(s[0], s0); s[1] = max(s[1], e0)
        zseq = {}
        for qn, (s0, e0) in spans.items():
            if e0 - s0 > 5_000_000:
                sys.stderr.write(f"      {qn}: lifted span too large, skipping validation\n")
                continue
            zseq[qn] = (s0, ucsc_seq(zbuild, qn, s0, e0))
        nmis = 0
        for g in guides:
            t = lifted.get(g["rank"])
            if not t:
                continue
            qn, s0, e0 = t
            g11 = ref[g["start"] - lo: g["start"] - lo + 23]
            if qn in zseq:
                base, s = zseq[qn]
                z = s[s0 - base: e0 - base]
                if z != g11 and rc(z) != g11:
                    g["zv9"] = "liftover-seq-mismatch"
                    nmis += 1
                    continue
            g["zv9"] = f"{qn}:{s0+1}-{e0}"
        nlift = sum(1 for g in guides if g["zv9"] and not g["zv9"].startswith("liftover"))
        sys.stderr.write(f"      {nlift}/{len(guides)} lifted and sequence-verified"
                         f"{f', {nmis} mismatched' if nmis else ''}\n")

        order = sorted(guides, key=lambda g: g["rank"])
        target = [g for g in (order if args.snp_top <= 0 else order[:args.snp_top])
                  if g["zv9"] and not g["zv9"].startswith("liftover")]
        wmap: dict[str, list[dict]] = {}
        for g in target:
            wmap.setdefault(g["zv9"], []).append(g)
        windows = list(wmap)
        sys.stderr.write(f"      querying SNPfisher for {len(windows)} windows "
                         f"(one batched request) ...\n")
        snp_fields = snpfisher_fields(cache)
        hits = snpfisher_batch(snp_fields, windows)
        for w, rows in hits.items():
            for g in wmap.get(w, []):
                g["snps"] = len(rows)
                g["snp_detail"] = rows
                for k, col in SNP_STRAINS.items():
                    g[f"snp_{k}"] = sum(1 for r in rows if _fq(r.get(col)) > 0)
        for g in target:
            if g["snps"] == "":
                g["snps"] = 0
                for k in SNP_STRAINS:
                    g[f"snp_{k}"] = 0
        nsnp = sum(1 for g in target if isinstance(g["snps"], int) and g["snps"] > 0)
        sys.stderr.write(f"      {nsnp} of {len(target)} carry at least one SNP\n")

    # 6. write -----------------------------------------------------------------
    sys.stderr.write("[6/6] writing outputs ...\n")
    exon_info = []
    if model:
        _t, _strand, _exons, _cds0, _cds1, _cdslen = model
        counts = {}
        for g in guides:
            if g.get("region") == "exon" and g.get("exon"):
                counts[g["exon"]] = counts.get(g["exon"], 0) + 1
        for i, (a, b) in enumerate(_exons, 1):
            coding = max(0, min(b, _cds1) - max(a, _cds0))
            exon_info.append({
                "n": i, "len": b - a, "coding": coding,
                "symmetry": ("noncoding" if coding == 0 else
                             "asymmetric" if coding % 3 else "symmetric"),
                "guides": counts.get(i, 0),
                "first": i == 1, "last": i == len(_exons)})
    stem = os.path.join(args.out, f"{gene}_{build}")
    write_csv(stem + ".csv", guides, job)
    write_html(stem + ".html", guides, job, gene, build, model, query, snp_fields,
               exon_info)
    sys.stderr.write(f"\n  CSV  {os.path.abspath(stem + '.csv')}\n")
    sys.stderr.write(f"  HTML {os.path.abspath(stem + '.html')}\n")
    sys.stderr.write(f"  job  {CHOPCHOP}/results/{job}/\n")

    if args.open:
        try:
            webbrowser.open("file://" + os.path.abspath(stem + ".html").replace(os.sep, "/"))
        except Exception as e:                      # never fail a good run over this
            print(f"  (could not open a browser: {e})", file=sys.stderr)


def write_csv(path: str, guides: list[dict], job: str) -> None:
    cols = ["rank", "target", "proto", "pam", "chrom", "start", "strand", "cut",
            "exon", "exon_total", "region", "symmetry", "bp_to_splice", "in_cds",
            "efficiency", "frameshift", "precision", "mh_del", "exp_indel",
            "gc", "selfcompl", "mm0", "mm1", "mm2", "mm3", "zv9", "snps", "snp_fli", "snp_tl", "snp_wik",
            "details_url"]
    nprim = max((len(g.get("primers") or []) for g in guides), default=0)
    pcols = []
    for i in range(1, min(nprim, 5) + 1):
        pcols += [f"primer{i}_fwd", f"primer{i}_rev", f"primer{i}_tm_f",
                  f"primer{i}_tm_r", f"primer{i}_size"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols + pcols, extrasaction="ignore")
        w.writeheader()
        for g in guides:
            row = dict(g)
            row["details_url"] = f"{CHOPCHOP}/results/{job}/details/{g['rank']}"
            for i, p in enumerate((g.get("primers") or [])[:5], 1):
                row[f"primer{i}_fwd"] = p["fwd"]; row[f"primer{i}_rev"] = p["rev"]
                row[f"primer{i}_tm_f"] = p["tm_f"]; row[f"primer{i}_tm_r"] = p["tm_r"]
                row[f"primer{i}_size"] = p["size"]
            w.writerow(row)


def write_html(path: str, guides, job, gene, build, model, query, snp_fields,
               exon_info=None) -> None:
    tpl_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "report_template.html")
    if not os.path.exists(tpl_path):
        raise PipelineError(f"missing HTML template at {tpl_path}",
                            "report_template.html must sit next to this script")
    tpl = open(tpl_path, encoding="utf-8").read()
    if "/*__DATA__*/" not in tpl:
        raise PipelineError("report_template.html has no /*__DATA__*/ placeholder")
    opts = query.get("opts") or []
    scoring = opts[opts.index("-scoringMethod") + 1] if "-scoringMethod" in opts else "unknown"
    payload = {
        "gene": gene, "build": build, "job": job,
        "jobUrl": f"{CHOPCHOP}/results/{job}/",
        "transcript": model[0] if model else None,
        "exonTotal": len(model[2]) if model else None,
        "cdsLen": model[5] if model else None,
        "aaLen": (model[5] // 3 - 1) if model else None,
        "scoring": scoring,
        "generated": time.strftime("%Y-%m-%d %H:%M"),
        "snpFields": snp_fields or {},
        "exonInfo": exon_info or [],
        "guides": guides,
    }
    blob = json.dumps(payload).replace("</", "<\\/")   # keep </script> out of the literal
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(tpl.replace("/*__DATA__*/", blob))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="CHOPCHOP -> annotated CSV + interactive HTML for zebrafish F0 knockouts.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--gene", help="gene symbol or Ensembl ID, e.g. yap1")
    ap.add_argument("--job-url", help="reuse an existing CHOPCHOP results URL instead of submitting")
    ap.add_argument("--build", default="danRer11", help="UCSC assembly (default danRer11)")
    ap.add_argument("--scoring", default="MORENO_MATEOS_2015",
                    help="CHOPCHOP efficiency model (default MORENO_MATEOS_2015, "
                         "the one Kroll recommends for zebrafish)")
    ap.add_argument("--pam", default="NGG")
    ap.add_argument("--product-min", type=int, default=150)
    ap.add_argument("--product-max", type=int, default=290)
    ap.add_argument("--snp-top", type=int, default=0,
                    help="query SNPfisher for the top N guides by CHOPCHOP rank. "
                         "Default 0 = all of them; SNPfisher is queried in a single "
                         "batched request, so restricting it saves almost nothing")
    ap.add_argument("--delay", type=float, default=0.2, help="seconds between per-guide fetches")
    ap.add_argument("--timeout", type=int, default=900, help="seconds to wait for the job")
    ap.add_argument("-o", "--out", default="chopchop_out", help="output directory")
    ap.add_argument("--open", action=argparse.BooleanOptionalAction, default=True,
                    help="open the HTML report in your default browser when done "
                         "(default: yes; use --no-open for scripts and batch runs)")
    args = ap.parse_args()
    if not args.gene and not args.job_url:
        ap.error("supply --gene or --job-url")
    try:
        run(args)
    except PipelineError as e:
        fail(str(e))
    except KeyboardInterrupt:
        fail("interrupted")


if __name__ == "__main__":
    main()
