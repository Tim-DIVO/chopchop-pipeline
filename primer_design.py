#!/usr/bin/env python3
"""
primer_design.py -- short, Tm-harmonised, SNP-free genotyping primers for a set
of CRISPR cut sites, sized for Illumina amplicon sequencing.

Why this exists
---------------
CHOPCHOP (and CRISPOR) hand you Primer3 pairs that are unusable for a MiSeq
amplicon panel, for three structural reasons:

  1. They design one locus at a time, so there is no shared annealing
     temperature across your guides. You cannot run the plate in one block.
  2. CRISPOR's `runPrimer3()` sets PRIMER_MIN/OPT/MAX_TM to tm-3 / tm / tm+3 and
     never sets PRIMER_PAIR_MAX_DIFF_TM, which Primer3 defaults to 100. Left and
     right primer can legally differ by 6 C.
  3. Neither excludes known variants from the primer footprint.

This script does the opposite: it treats the whole guide set as one design
problem, scans for the single Tm that satisfies the most loci, and only then
relaxes anything -- reporting exactly what it relaxed, per locus.

Deterministic. No judgement calls live in here; the interpretation layer is the
`crispr-primers` skill (see SKILL.md). Everything this script decides is a
documented rule with a flag attached to it.

Design targets (provenance in README / SKILL.md):
  * product <= 233 bp     -- the user's MiSeq budget; also ~67 bp of read overlap
                             on a 2x150 run, comfortably over Illumina's >=50 bp
  * cut >= 40 bp from either product end
                          -- matches Kroll et al. 2021's own MiSeq panel
                             (42 loci: cut 39-93 bp from the nearer end), and
                             clears CRISPResso2's default 15 bp end-exclusion
  * primers never overlap the protospacer+PAM
  * left/right Tm within 0.5 C, and one Tm for the entire plate

Usage
-----
  python primer_design.py \
      --from out_geneA/geneA_danRer11.csv:3,6,34,15,33,5 \
      --from out/yap1_danRer11.csv:5,16,17,13,10,21 \
      -o primers_out

  # already know the Tm you want to run at:
  python primer_design.py --from ... --tm 60.0

Outputs <out>/primers.csv, <out>/order_sheet.csv and <out>/primers.html.
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

# --------------------------------------------------------------------------- SSL
# The conda cert store on this machine is missing UCSC's issuer; without this,
# every api.genome.ucsc.edu call dies with CERTIFICATE_VERIFY_FAILED.
try:
    import ssl

    import certifi

    _SSL_CTX = ssl.create_default_context(cafile=certifi.where())
    ssl._create_default_https_context = lambda *a, **k: _SSL_CTX  # noqa: SLF001
except Exception:  # pragma: no cover - fall back to the system store
    _SSL_CTX = None

try:
    import primer3
except ImportError:
    sys.exit("primer3-py is required:  pip install primer3-py")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chopchop_pipeline import (  # noqa: E402
    ensure_chain,
    lift,
    load_chain,
    rc,
    snpfisher_batch,
    snpfisher_fields,
    ucsc_seq,
)

UA = "Mozilla/5.0 (chopchop-pipeline primer_design)"
ENSEMBL = "https://rest.ensembl.org"
UCSC_PCR = "https://genome.ucsc.edu/cgi-bin/hgPcr"

# Illumina 16S Metagenomic Sequencing Library Prep, doc 15044223 Rev. B, p.3.
# 33 nt + 34 nt = 67. These are NOT part of the 233 bp budget: the sequencing
# read starts after the overhang, so they only lengthen the flowcell fragment.
OVERHANG_F = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
OVERHANG_R = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"


# ========================================================================== util
def log(msg: str) -> None:
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def get_json(url: str, timeout: int = 60):
    req = urllib.request.Request(url, headers={"User-Agent": UA,
                                               "Accept": "application/json"})
    return json.loads(urllib.request.urlopen(req, timeout=timeout).read().decode())


def get_text(url: str, timeout: int = 120) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    return urllib.request.urlopen(req, timeout=timeout).read().decode("utf-8", "replace")


# ======================================================================== inputs
class Locus:
    """One cut site to design around."""

    def __init__(self, gene: str, rank: str, chrom: str, cut: int,
                 target: str, strand: str, exon: str = "", aa: str = ""):
        self.gene, self.rank = gene, rank
        self.chrom, self.cut = chrom, cut
        self.target, self.strand = target.upper(), strand
        self.exon, self.aa = exon, aa
        self.name = f"{gene}_g{rank}"
        # filled in later
        self.template = ""
        self.tpl_start = 0          # 0-based genomic coord of template[0]
        self.guide_lo = self.guide_hi = 0   # protospacer+PAM span, template coords
        self.cut_off = 0            # cut position in template coords
        self.snps: dict[int, list] = {}     # template offset -> variant records
        self.pairs: list[dict] = []
        self.tier = None
        self.notes: list[str] = []
        #: Filled once, after design, by recommend_pairs -- and then used by
        #: every writer. Recomputing it per writer let two call sites disagree
        #: on a borderline overlap, because they measured spans in different
        #: coordinate systems.
        self.recommended: list[int] = []
        self.rec_why: str = ""
        self.rec_independent: bool = False
        #: False once a variant lookup has failed. A locus with snp_ok False is
        #: reported as UNVERIFIED no matter which tier it reached, because the
        #: masking it was designed against was incomplete.
        self.snp_ok = True


def read_loci(spec: str) -> list[Locus]:
    """--from path/to.csv:3,6,34   (ranks are CHOPCHOP ranks in that CSV)"""
    if ":" not in spec:
        sys.exit(f"--from needs the form FILE.csv:rank,rank,...  got {spec!r}")
    path, ranks = spec.rsplit(":", 1)
    wanted = [r.strip() for r in ranks.split(",") if r.strip()]
    if not os.path.exists(path):
        sys.exit(f"no such file: {path}")
    gene = os.path.basename(path).split("_")[0]
    rows = {r["rank"]: r for r in csv.DictReader(open(path, newline=""))}
    out = []
    for rk in wanted:
        r = rows.get(rk)
        if r is None:
            sys.exit(f"rank {rk} not found in {path}")
        out.append(Locus(gene, rk, r["chrom"], int(r["cut"]), r["target"],
                         r["strand"], r.get("exon", ""), r.get("aa", "")))
    return out


# =================================================================== template
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".cache", "seq")


def cached_seq(build: str, chrom: str, lo: int, hi: int, retries: int = 5) -> str:
    """UCSC sequence, cached on disk and backed off politely.

    Without the cache, every re-run re-fetches the same windows and UCSC starts
    answering 429. The sequence for a fixed assembly and interval never changes,
    so there is no reason to ask twice.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, f"{build}_{chrom}_{lo}_{hi}.txt")
    if os.path.exists(path):
        s = open(path).read().strip()
        if s:
            return s
    wait = 5.0
    for attempt in range(retries):
        try:
            s = ucsc_seq(build, chrom, lo, hi)
            with open(path, "w") as fh:
                fh.write(s)
            return s
        except Exception as e:
            if "429" not in str(e) or attempt == retries - 1:
                raise
            log(f"      UCSC rate-limited; waiting {wait:.0f}s "
                f"(attempt {attempt + 1}/{retries})")
            time.sleep(wait)
            wait *= 2
    raise RuntimeError("unreachable")


def fetch_template(loc: Locus, build: str, flank: int) -> None:
    lo = max(0, loc.cut - flank)
    hi = loc.cut + flank
    loc.template = cached_seq(build, loc.chrom, lo, hi)
    loc.tpl_start = lo

    # Anchor on the 23-mer itself rather than trusting coordinate conventions.
    t = loc.target
    i = loc.template.find(t)
    if i < 0:
        i = loc.template.find(rc(t))
    if i < 0:
        raise RuntimeError(f"{loc.name}: protospacer+PAM not found in the fetched "
                           f"window -- check build/coordinates")
    loc.guide_lo, loc.guide_hi = i, i + len(t)
    loc.cut_off = loc.cut - lo


# ====================================================================== variants
def ensembl_variants(loc: Locus, retries: int = 5) -> list[dict]:
    """dbSNP/Ensembl variation on GRCz11 (== danRer11). 1-based, no 'chr'.

    Cached and retried. On final failure this RAISES rather than returning an
    empty list: "the variant service was down" and "there are no variants here"
    must never produce the same primers. The caller marks the locus unverified.
    """
    chrom = loc.chrom.replace("chr", "")
    lo1 = loc.tpl_start + 1
    hi1 = loc.tpl_start + len(loc.template)
    os.makedirs(CACHE_DIR, exist_ok=True)
    cpath = os.path.join(CACHE_DIR, f"ens_{chrom}_{lo1}_{hi1}.json")
    if os.path.exists(cpath):
        try:
            recs = json.load(open(cpath))
        except Exception:
            recs = None
        if recs is not None:
            return _ens_offsets(recs, loc)

    url = (f"{ENSEMBL}/overlap/region/danio_rerio/{chrom}:{lo1}-{hi1}"
           f"?feature=variation;content-type=application/json")
    wait = 5.0
    recs = None
    for attempt in range(retries):
        try:
            recs = get_json(url, timeout=90)
            break
        except Exception as e:
            if attempt == retries - 1:
                raise RuntimeError(
                    f"Ensembl variation lookup failed for {loc.name} after "
                    f"{retries} attempts: {type(e).__name__}") from e
            log(f"      Ensembl busy ({type(e).__name__}); waiting {wait:.0f}s "
                f"(attempt {attempt + 1}/{retries})")
            time.sleep(wait)
            wait *= 2
    with open(cpath, "w") as fh:
        json.dump(recs, fh)
    return _ens_offsets(recs, loc)


def _ens_offsets(recs, loc: Locus) -> list[dict]:
    out = []
    for v in recs:
        off = int(v["start"]) - 1 - loc.tpl_start
        if 0 <= off < len(loc.template):
            out.append(dict(off=off, id=v.get("id", ""), source="ensembl",
                            alleles="/".join(v.get("alleles", []))))
    return out


def snpfisher_variants(loci: list[Locus], build: str, cache: str) -> None:
    """Strain SNPs (TL/WIK/FLI) -- Zv9 only, so each window has to be lifted."""
    chain = ensure_chain(build, cache)
    if not chain:
        log("      ! no liftOver chain for this build; skipping SNPfisher")
        return
    try:
        fields = snpfisher_fields(cache)
    except Exception as e:
        log(f"      ! SNPfisher form unavailable ({type(e).__name__}); skipping")
        return

    windows, back = [], {}
    by_chrom: dict[str, list[Locus]] = {}
    for lo in loci:
        by_chrom.setdefault(lo.chrom, []).append(lo)

    for chrom, group in by_chrom.items():
        blocks = load_chain(chain, chrom)
        for loc in group:
            # Lift every base of the window individually and keep the exact
            # inverse map. Interpolating across the window would be wrong
            # wherever the chain has an indel, which is the whole point of a
            # chain file.
            inv: dict[int, int] = {}
            zchrom = None
            for off in range(len(loc.template)):
                hit = lift(blocks, loc.tpl_start + off)
                if not hit:
                    continue
                if zchrom is None:
                    zchrom = hit[0]
                elif hit[0] != zchrom:
                    continue
                inv[hit[1]] = off
            if not inv:
                loc.notes.append("no Zv9 liftOver; SNPfisher not checked")
                continue
            z0, z1 = min(inv), max(inv)
            # SNPfisher windows and its Pos column are 1-based (the pipeline
            # builds them as s0+1 at chopchop_pipeline.py:806); `inv` is keyed
            # on the 0-based lifted coordinate.
            key = f"{zchrom}:{z0 + 1}-{z1 + 1}"
            windows.append(key)
            back[key] = (loc, inv)
    if not windows:
        return
    try:
        res = snpfisher_batch(fields, windows)
    except Exception as e:
        log(f"      ! SNPfisher query failed ({type(e).__name__}); skipping")
        return
    for key, rows in res.items():
        if key not in back:
            continue
        loc, inv = back[key]
        for r in rows:
            try:
                zp = int(r["Pos"])
            except (KeyError, ValueError):
                continue
            off = inv.get(zp - 1)
            if off is not None:
                loc.snps.setdefault(off, []).append(
                    dict(off=off, id=f"{r.get('Chr','')}:{zp}", source="snpfisher",
                         alleles=f"{r.get('Ref','')}/{r.get('Mut','')}",
                         fli=r.get("FLI_FQ", ""), tl=r.get("TL_FQ", ""),
                         wik=r.get("WIK_FQ", "")))


# ======================================================================= primer3
#: The relaxation ladder. Earlier tiers are strictly preferable. Tiers 1-3 keep
#: the shared annealing temperature intact and only give up on other things;
#: only tiers 4-5 break Tm uniformity, because that is what costs you the
#: single-block plate.
#: `pmin_drop` and `flank_drop` are subtracted from the user's --min-product /
#: --min-flank, so the ladder keeps its shape whatever those are set to.
TIERS = [
    dict(name="strict", tm_halfwin=0.5, diff_tm=0.5, pmin_drop=0, flank_drop=0,
         snp_mode="all",
         why="all constraints met"),
    dict(name="snp-3prime", tm_halfwin=0.5, diff_tm=0.5, pmin_drop=0, flank_drop=0,
         snp_mode="three_prime",
         why="a SNP is tolerated in the primer 5' region, but never in the 3' 8 nt"),
    dict(name="flank-30", tm_halfwin=0.5, diff_tm=0.5, pmin_drop=0, flank_drop=10,
         snp_mode="three_prime",
         why="cut allowed 10 bp closer to a product end"),
    dict(name="product-100", tm_halfwin=0.5, diff_tm=0.5, pmin_drop=40, flank_drop=10,
         snp_mode="three_prime",
         why="shorter products allowed"),
    dict(name="tm-1.0", tm_halfwin=1.0, diff_tm=0.75, pmin_drop=40, flank_drop=10,
         snp_mode="three_prime",
         why="Tm window widened to +/-1.0 C -- still one block, but less uniform"),
    dict(name="tm-1.5", tm_halfwin=1.5, diff_tm=1.0, pmin_drop=40, flank_drop=15,
         snp_mode="three_prime",
         why="Tm window widened to +/-1.5 C -- check this locus on a gradient"),
]

THREE_PRIME_NT = 8   # a SNP this close to the 3' end kills extension

# --------------------------------------------------------------- length policy
# Primer length is the cheapest specificity available here. An exact 20-mer is
# expected to occur by chance about 4**5 = 1024 times more often in a 1.4 Gb
# genome than an exact 25-mer, and lengthening an oligo costs nothing but a few
# cents of synthesis. So push length as far as the other constraints allow.
#
# Be honest about the limit of that argument: the Tm window is fixed (one block
# for the whole plate), so a longer primer is necessarily more AT-rich, and
# AT-rich primers tolerate a *mismatched* partial site about as well as the
# short GC-rich ones do. Length buys uniqueness of the exact match, not
# thermodynamic discrimination -- which is why it is a default, not a
# substitute for blast_specificity.py.
#
# The defaults are the measured optimum on a twelve-guide test set, not a guess.
# Sweeping (min, opt, max) over the twelve loci (scratch experiment, 2026-11-16):
#
#   18/20/27  Tm 57.75  mean oligo 20.6 nt  11 strict + 1 flank-30   6/12 independent
#   18/23/27  Tm 59.75  mean oligo 22.9 nt  11 strict + 1 flank-30   7/12 independent
#   18/25/27  Tm 59.75  mean oligo 23.0 nt  11 strict + 1 flank-30   5/12 independent
#   20/25/27  Tm 59.75  mean oligo 23.2 nt  10 strict, hardest locus -> snp-3prime
#   22/26/30  Tm 61.50  mean oligo 24.6 nt   9 strict, hardest locus FAILED entirely
#
# Past opt 23 the mean length stops improving (23.0 vs 22.9) while the candidate
# pool thins, because a longer primer inside a fixed Tm window must be AT-rich
# and there are fewer places to put it -- which costs independent backups and,
# past min 20, costs whole loci. So 23, and raising `min` above 18 is a trap:
# the hardest locus in the test set needs one short GC-rich oligo and
# fails without it.
SIZE = dict(min=18, opt=23, max=27)

# How the shortlist is ordered. "length" = specificity-first (see pair_sort_key);
# "penalty" = the old behaviour, Primer3's own idea of a good pair.
SORT_MODE = "length"

#: Product-size optimum, in bp of genomic amplicon (before the MiSeq overhangs).
#: Set from --opt-product. The MiSeq amplicon protocol used on the fish floor
#: works best around 200 bp, which becomes ~267 bp once the two Illumina tails
#: are added; much shorter than that wastes read length and much longer loses
#: cluster density. This is a *preference*, not a constraint: --min-product and
#: --max-product are the hard bounds.
PRODUCT_OPT = 200

# Set from --specificity so write_candidates() can record how specificity was
# measured; None means it was not.
SPEC_META = None


def pair_sort_key(p: dict):
    """Rank one locus's candidate pairs.

    Two quantities are deliberately bucketed, for the same reason: a difference
    far inside the noise of what it predicts must not outrank a difference that
    matters.

    - dTm to 0.1 C. Every pair here already sits inside the tier's dTm ceiling,
      and on a single annealing block a pair differing by 0.03 C is not a better
      experiment than one differing by 0.08 C. Sorting on that second decimal
      let noise outrank primer length.
    - the shorter oligo's length to 2 nt. One base changes the chance of a
      spurious exact match by 4x; a 59 bp longer amplicon is worth more than
      that, so a 1 nt edge should not win on its own.
    - distance from the product optimum to 10 bp. 196 bp and 203 bp are the
      same experiment on a MiSeq; ranking on single bp there would let product
      size outrank primer length, which is the quantity that actually predicts
      a spurious match.

    Order: **predicted specificity first** when it is known (a pair that
    amplifies a second locus does not genotype anything, however pretty its Tm),
    then dTm bucket, then **variants under either primer** (a SNP there risks
    allele dropout in exactly the genotyping assay these primers exist for, so
    it outranks every length argument), then the shorter oligo's length bucket,
    then **closeness to the product optimum**, then the finer length detail,
    then Primer3.

    A pair whose specificity was never measured sorts between "clean" and
    "has an off-target": not promoted over a pair known to be clean, not
    condemned on an absence of evidence.
    """
    if SORT_MODE == "penalty":
        return (p["dtm"], p["penalty"])
    lo = min(len(p["fwd"]), len(p["rev"]))
    off = p.get("ispcr_off")
    spec = 0 if off == 0 else (1 if off is None else 1 + min(off, 50))
    return (spec,
            round(p["dtm"] / 0.1),
            p["snp_fwd"] + p["snp_rev"],
            -(lo // 2),
            abs(p["size"] - PRODUCT_OPT) // 10,
            -lo,
            -(len(p["fwd"]) + len(p["rev"])),
            abs(p["size"] - PRODUCT_OPT),
            p["penalty"])


def _excluded(loc: Locus, mode: str) -> list[list[int]]:
    """SEQUENCE_EXCLUDED_REGION: primers may not overlap these."""
    if mode == "none":
        return []
    return [[off, 1] for off in sorted(loc.snps)]


def design_locus(loc: Locus, tm: float, tier: dict, conc: dict,
                 num_return: int = 20) -> list[dict]:
    half, flank = tier["tm_halfwin"], tier["flank"]

    # SEQUENCE_TARGET forces the product to contain this span. Setting it to the
    # guide plus `flank` on both sides is what guarantees the cut sits at least
    # `flank` bp inside each product end, and that no primer lands on the guide.
    lo = max(0, min(loc.guide_lo, loc.cut_off - flank))
    hi = min(len(loc.template), max(loc.guide_hi, loc.cut_off + flank))

    seq_args = {
        "SEQUENCE_ID": loc.name,
        "SEQUENCE_TEMPLATE": loc.template,
        "SEQUENCE_TARGET": [lo, hi - lo],
    }
    if tier["snp_mode"] == "all":
        excl = _excluded(loc, "all")
        if excl:
            seq_args["SEQUENCE_EXCLUDED_REGION"] = excl

    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_PRODUCT_SIZE_RANGE": [[tier["pmin"], tier["pmax"]]],
        "PRIMER_PRODUCT_OPT_SIZE": tier["popt"],
        # Primer3 ignores PRODUCT_OPT_SIZE unless these weights are non-zero,
        # and both default to 0. Undershooting is weighted 5x harder than
        # overshooting: a product below ~175 bp starts to crowd the cut site
        # against the read ends, while a long one only costs efficiency.
        "PRIMER_PAIR_WT_PRODUCT_SIZE_LT": 0.5,
        "PRIMER_PAIR_WT_PRODUCT_SIZE_GT": 0.1,
        "PRIMER_MIN_TM": tm - half,
        "PRIMER_OPT_TM": tm,
        "PRIMER_MAX_TM": tm + half,
        "PRIMER_PAIR_MAX_DIFF_TM": tier["diff_tm"],
        "PRIMER_MIN_SIZE": SIZE["min"],
        "PRIMER_OPT_SIZE": SIZE["opt"],
        "PRIMER_MAX_SIZE": SIZE["max"],
        "PRIMER_MIN_GC": 35.0,
        "PRIMER_MAX_GC": 65.0,
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_MAX_SELF_ANY_TH": 45.0,
        "PRIMER_MAX_SELF_END_TH": 35.0,
        "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45.0,
        "PRIMER_PAIR_MAX_COMPL_END_TH": 35.0,
        "PRIMER_MAX_HAIRPIN_TH": 40.0,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_GC_CLAMP": 1,
        "PRIMER_SALT_MONOVALENT": conc["mv"],
        "PRIMER_SALT_DIVALENT": conc["dv"],
        "PRIMER_DNTP_CONC": conc["dntp"],
        "PRIMER_DNA_CONC": conc["dna"],
    }

    try:
        res = primer3.design_primers(seq_args, global_args)
    except Exception as e:
        loc.notes.append(f"primer3 error: {type(e).__name__}: {e}")
        return []

    n = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)
    out = []
    for i in range(n):
        f_seq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
        r_seq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        f_start, f_len = res[f"PRIMER_LEFT_{i}"]
        r_start, r_len = res[f"PRIMER_RIGHT_{i}"]
        prod_lo = f_start
        prod_hi = r_start + 1
        size = int(res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"])

        # SNPs under each primer, split by how damaging the position is
        f_span = range(f_start, f_start + f_len)
        r_span = range(r_start - r_len + 1, r_start + 1)
        f_snp = [o for o in loc.snps if o in f_span]
        r_snp = [o for o in loc.snps if o in r_span]
        # 3' end of the left primer is its highest coordinate; of the right, its lowest
        f_3p = [o for o in f_snp if o >= f_start + f_len - THREE_PRIME_NT]
        r_3p = [o for o in r_snp if o <= r_start - r_len + THREE_PRIME_NT]

        if tier["snp_mode"] == "three_prime" and (f_3p or r_3p):
            continue   # never acceptable

        pair = dict(
            fwd=f_seq, rev=r_seq,
            tm_f=round(float(res[f"PRIMER_LEFT_{i}_TM"]), 2),
            tm_r=round(float(res[f"PRIMER_RIGHT_{i}_TM"]), 2),
            gc_f=round(float(res[f"PRIMER_LEFT_{i}_GC_PERCENT"]), 1),
            gc_r=round(float(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"]), 1),
            size=size,
            cut_from_left=loc.cut_off - prod_lo,
            cut_from_right=prod_hi - loc.cut_off,
            snp_fwd=len(f_snp), snp_rev=len(r_snp),
            snp_3prime=len(f_3p) + len(r_3p),
            penalty=round(float(res[f"PRIMER_PAIR_{i}_PENALTY"]), 3),
            prod_lo=prod_lo, prod_hi=prod_hi,
        )
        pair["dtm"] = round(abs(pair["tm_f"] - pair["tm_r"]), 2)
        try:
            pair["hetero_dg"] = round(
                primer3.calc_heterodimer(f_seq, r_seq, mv_conc=conc["mv"],
                                         dv_conc=conc["dv"], dntp_conc=conc["dntp"],
                                         dna_conc=conc["dna"]).dg / 1000.0, 2)
        except Exception:
            pair["hetero_dg"] = None
        out.append(pair)
    return out


# ================================================================== specificity
def ispcr_local(ispcr_bin: str, genome: str, pairs: list[tuple[str, str, str]],
                max_size: int = 4000) -> dict[str, list[str]]:
    """Genome-wide in-silico PCR with UCSC's `isPcr` binary against a local
    .2bit/.fa. This is the sanctioned way to do this at scale.

    NOT done by scraping genome.ucsc.edu/cgi-bin/hgPcr: that endpoint sits
    behind a Cloudflare Turnstile bot check which explicitly points automated
    users at https://genome.ucsc.edu/FAQ/FAQdownloads.html#CAPTCHA. Getting
    past it means lying about the User-Agent, so the tool does not do it.
    """
    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        pf = os.path.join(td, "pairs.txt")
        with open(pf, "w") as fh:
            for name, f, r in pairs:
                fh.write(f"{name}\t{f}\t{r}\n")
        out = os.path.join(td, "out.fa")
        cmd = [ispcr_bin, genome, pf, f"-maxSize={max_size}", "-out=fa", out]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=1800)
        except Exception as e:
            log(f"      ! local isPcr failed: {type(e).__name__}: {e}")
            return {}
        hits: dict[str, list[str]] = {name: [] for name, _, _ in pairs}
        for line in open(out):
            if not line.startswith(">"):
                continue
            # >chr20:36408662+36408837 geneA_g3 176bp FWD REV
            parts = line[1:].split()
            if len(parts) >= 3:
                hits.setdefault(parts[1], []).append(f"{parts[0]} {parts[2]}")
        return hits


def _ispcr_web_disabled(build: str, fwd: str, rev: str, max_size: int = 4000,
                        org: str = "Zebrafish") -> tuple[int, list[str]]:
    """Kept only to document why there is no web fallback. Do not re-enable."""
    q = dict(org=org, db=build, wp_target="genome", wp_f=fwd, wp_r=rev,
             wp_size=str(max_size), wp_perfect="15", wp_good="15", Submit="submit")
    try:
        txt = get_text(f"{UCSC_PCR}?{urllib.parse.urlencode(q)}", timeout=180)
    except Exception as e:
        return -1, [f"isPCR failed: {type(e).__name__}"]
    # Results live inside a <PRE> block as FASTA; the header line of each record
    # is one amplicon: ">chr20:36408662+36408837 176bp FWD REV". Scanning the raw
    # HTML for lines starting with ">" silently finds nothing, because the first
    # one is glued to the opening <PRE> tag.
    hits = []
    for block in re.findall(r"<PRE>(.*?)</PRE>", txt, re.S | re.I):
        for line in re.sub(r"<[^>]+>", "", block).splitlines():
            line = line.strip()
            if not line.startswith(">"):
                continue
            m = re.match(r">(\S+)\s+(\d+bp)", line)
            hits.append(f"{m.group(1)} {m.group(2)}" if m else line[1:81])
    return len(hits), hits


# ===================================================================== Tm scan
def scan_tm(loci: list[Locus], lo: float, hi: float, step: float,
            conc: dict) -> list[dict]:
    """How many loci get a strict-tier pair at each candidate Tm?"""
    tier = TIERS[0]
    rows = []
    tm = lo
    while tm <= hi + 1e-9:
        ok, pens = 0, []
        for loc in loci:
            # 5 is deliberate: the scan only asks *whether* a locus solves at
            # this Tm, not which pair wins, and it runs this for every step of
            # the sweep across every locus.
            pairs = design_locus(loc, tm, tier, conc, num_return=5)
            if pairs:
                ok += 1
                pens.append(min(p["penalty"] for p in pairs))
        rows.append(dict(tm=round(tm, 2), n_ok=ok,
                         mean_penalty=round(sum(pens) / len(pens), 3) if pens else None))
        tm += step
    return rows


def pick_tm(scan: list[dict]) -> float:
    best = max(scan, key=lambda r: (r["n_ok"],
                                    -(r["mean_penalty"] if r["mean_penalty"] is not None else 1e9)))
    return best["tm"]


# ======================================================================= output
def recommend_pairs(pairs: list[dict]) -> tuple[list[int], str, bool]:
    """Pick a primary pair and a backup.

    The primary is rank 1 -- the list is already sorted by dTm then Primer3
    penalty, so that is the best pair on the criteria the search optimised.

    The backup is the interesting choice. Primer3 typically holds one side fixed
    and walks the other, so consecutive ranks often re-use a primer. A backup
    sharing a primer with the primary inherits its failure: if the reverse
    primer is the reason the PCR did not work, the "backup" fails identically.
    So prefer a pair that shares neither side.

    Two pairs per guide are always recommended when two exist, because that is
    the plate layout. Where no independent pair exists the near-duplicate is
    still returned, with `independent=False` and a `why` that says exactly what
    is and is not insured -- rather than silently handing back one pair and
    leaving a hole in the order.

    Among independent candidates, rank by: no variant anywhere under a primer,
    then no strong self-dimer, then dTm, then Primer3 penalty.
    """
    if not pairs:
        return [], "no pairs at any tier", False
    first = pairs[0]
    rest = pairs[1:]
    if not rest:
        return [first["i"]], "only one pair available -- no backup exists", False

    clean = lambda p: p["snp_3prime"] == 0          # noqa: E731

    def spans(p):
        """Coordinates each primer actually binds.

        Called both on raw Primer3 output (template offsets, prod_lo/prod_hi)
        and on the dicts written into candidates.json (genomic, f_start/f_end).
        Either is fine -- only relative overlap within one locus is used.
        """
        if "prod_lo" in p:
            return ((p["prod_lo"], p["prod_lo"] + len(p["fwd"])),
                    (p["prod_hi"] - len(p["rev"]), p["prod_hi"]))
        return ((p["f_start"], p["f_end"]), (p["r_start"], p["r_end"]))

    def overlap(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))

    def independent(p):
        """Do neither of this pair's primers sit on the primary's sites?

        Comparing sequences is not enough: Primer3 routinely returns a primer
        that is its neighbour extended or trimmed by a base or two. Those are
        different strings but the same binding site, and they share every
        failure mode that matters -- a SNP under them, local secondary
        structure, a mispriming partner. Compare coordinates instead, and call
        it the same site when they overlap by more than half the shorter oligo.
        """
        fa, ra = spans(first)
        fb, rb = spans(p)
        f_ov = overlap(fa, fb) / min(fa[1] - fa[0], fb[1] - fb[0])
        r_ov = overlap(ra, rb) / min(ra[1] - ra[0], rb[1] - rb[0])
        return f_ov <= 0.5 and r_ov <= 0.5

    indep = [p for p in rest if independent(p) and clean(p)]

    def key(p):
        dg = p.get("hetero_dg")
        off = p.get("ispcr_off")
        spec = 0 if off == 0 else (1 if off is None else 1 + min(off, 50))
        return (spec,
                p["snp_fwd"] + p["snp_rev"],
                0 if (dg is None or dg > -9.0) else 1,
                p["dtm"], p["penalty"])

    if indep:
        b = min(indep, key=key)
        why = (f"#{b['i']} shares neither primer with #{first['i']}, so it is a "
               f"genuine independent backup")
        independent = True
    else:
        # No candidate binds a different site, so the best available backup is a
        # near-duplicate: the same primer shifted a base or two. It shares the
        # failure modes it is meant to insure against (a variant under the site,
        # local structure, a mispriming partner), so it is weak insurance -- but
        # two pairs per guide is the plate layout, so order it anyway and label
        # it. A backup that fails for the same reason is still a second shot at
        # a primer that was simply synthesised badly or pipetted wrong.
        b = min([p for p in rest if clean(p)] or rest, key=key)
        fa, ra = spans(first)
        fb, rb = spans(b)
        shared = []
        if overlap(fa, fb) / min(fa[1] - fa[0], fb[1] - fb[0]) > 0.5:
            shared.append("forward")
        if overlap(ra, rb) / min(ra[1] - ra[0], rb[1] - rb[0]) > 0.5:
            shared.append("reverse")
        why = (f"#{b['i']} is the best backup available, but it re-uses the "
               f"{' and '.join(shared) or 'same'} primer site, so it shares that "
               f"site's failure modes -- weak insurance, not independent. No "
               f"candidate at these settings binds a different site; widen the "
               f"search (longer product, smaller min-flank) if a genuinely "
               f"independent backup matters here")
        return [first["i"], b["i"]], why, False
    # Returned as a flag, not left to be inferred from the prose: "no
    # independent pair exists" contains the word "independent", and a caller
    # string-matching on it reads the exact opposite of the truth.
    return [first["i"], b["i"]], why, independent


def write_candidates(path: str, loci: list[Locus], tm: float, top: int,
                     args=None, scan=None) -> None:
    """Dump the ranked candidate pairs per locus, for the HTML report.

    `write_csv` emits only `pairs[0]`, the winner. That is the right thing for an
    order sheet but it hides the decision: a reader cannot see what the runner-up
    looked like or how close the call was. This writes the whole shortlist with
    genomic coordinates, so the report can lay the options out the way CHOPCHOP's
    detail page does and the choice can be argued with.

    Coordinates are emitted 1-based inclusive to match what CHOPCHOP displays,
    even though everything internal here is 0-based -- the conversion happens
    once, here, rather than in the JavaScript.
    """
    out = []
    for loc in loci:
        g1 = loc.tpl_start + loc.guide_lo + 1
        rec = {
            "name": loc.name, "gene": loc.gene, "rank": loc.rank,
            "chrom": loc.chrom, "strand": loc.strand, "cut": loc.cut,
            "exon": loc.exon, "aa": loc.aa,
            "tier": (loc.tier["name"] if loc.tier else "FAILED") if loc.snp_ok
                    else "UNVERIFIED-" + (loc.tier["name"] if loc.tier else "FAILED"),
            "tier_why": (loc.tier["why"] if loc.tier else "no pair at any tier"),
            "snp_ok": loc.snp_ok,
            "shared_tm": tm,
            "target": loc.target,
            "guide_start": g1,
            "guide_end": loc.tpl_start + loc.guide_hi,
            "tpl_start": loc.tpl_start + 1,
            "template": loc.template,
            # variant positions as 1-based genomic coords, for the track
            "snps": sorted(loc.tpl_start + o + 1 for o in loc.snps),
            "notes": list(loc.notes),
            "pairs": [],
        }
        keep = list(loc.pairs[:top])
        have = {p["i"] for p in keep}
        for i in loc.recommended:
            if i not in have:
                extra = next((p for p in loc.pairs if p["i"] == i), None)
                if extra:
                    keep.append(extra)
        keep.sort(key=lambda p: p["i"])
        for p in keep:
            f_lo = loc.tpl_start + p["prod_lo"]
            f_hi = f_lo + len(p["fwd"])
            r_hi = loc.tpl_start + p["prod_hi"]
            r_lo = r_hi - len(p["rev"])
            rec["pairs"].append({
                # the rank assigned at sort time -- NOT the position in this
                # list, which can differ once a deep recommended pair is
                # appended to the shortlist
                "i": p["i"],
                "fwd": p["fwd"], "rev": p["rev"],
                "tm_f": p["tm_f"], "tm_r": p["tm_r"], "dtm": p["dtm"],
                "gc_f": p["gc_f"], "gc_r": p["gc_r"],
                "size": p["size"],
                "cut_from_left": p["cut_from_left"],
                "cut_from_right": p["cut_from_right"],
                "snp_fwd": p["snp_fwd"], "snp_rev": p["snp_rev"],
                "snp_3prime": p["snp_3prime"],
                "penalty": p["penalty"], "hetero_dg": p.get("hetero_dg"),
                "ispcr_off": p.get("ispcr_off"),
                "ispcr_verdict": p.get("ispcr_verdict"),
                "ispcr_hits": p.get("ispcr_hits") or [],
                "f_start": f_lo + 1, "f_end": f_hi,
                "r_start": r_lo + 1, "r_end": r_hi,
                "prod_start": f_lo + 1, "prod_end": r_hi,
                "ispcr_n": p.get("ispcr_n"),
            })
        rec["recommended"] = list(loc.recommended)
        rec["rec_why"] = loc.rec_why
        rec["rec_independent"] = loc.rec_independent
        # Fingerprint of this locus's shortlist and recommendation. The report
        # stores hand-ticked selections in localStorage by rank, and ranks mean
        # nothing across a re-run: if the design changed, rank 2 is a different
        # pair of oligos. Stale picks used to survive and silently mask the new
        # recommendation. The report keeps this signature next to the stored
        # picks and drops them when it no longer matches.
        rec["sel_sig"] = hashlib.sha1(json.dumps(
            [[p["i"], p["fwd"], p["rev"]] for p in rec["pairs"]]
            + [list(loc.recommended)], sort_keys=True).encode()
        ).hexdigest()[:12]
        out.append(rec)
    # Record the conditions the search ran under. Without them the report shows
    # numbers with no way to tell what produced them, and no way to re-run with
    # one thing changed -- which is the whole point of the panel that reads this.
    params = None
    if args is not None:
        params = {
            "sources": list(args.sources),
            "build": args.build,
            "out": args.out,
            "n_loci": len(loci),   # the scan is over ALL loci, not one gene
            "tm_mode": "fixed" if args.tm is not None else "scan",
            "tm_chosen": tm,
            "tm_lo": args.tm_lo, "tm_hi": args.tm_hi, "tm_step": args.tm_step,
            "max_product": args.max_product, "min_product": args.min_product,
            "min_flank": args.min_flank, "flank_window": args.flank_window,
            "snpfisher": not args.no_snpfisher,
            "max_candidates": top,
            "num_return": args.num_return,
            "mv": args.mv, "dv": args.dv, "dntp": args.dntp, "dna": args.dna,
            "three_prime_nt": THREE_PRIME_NT,
            "primer_size": [SIZE["min"], SIZE["opt"], SIZE["max"]],
            "sort": SORT_MODE,
            "opt_product": args.opt_product,
            "gc_range": [35.0, 65.0],
            "diff_tm_strict": TIERS[0]["diff_tm"],
            "tiers": [{"name": t["name"], "why": t["why"]} for t in TIERS],
            "overhang_f": OVERHANG_F, "overhang_r": OVERHANG_R,
            "generated": time.strftime("%Y-%m-%d %H:%M"),
            "primer3": _primer3_version(),
            "specificity": SPEC_META,
        }
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"shared_tm": tm, "params": params,
                   "tm_scan": scan or [], "loci": out}, fh)


def _primer3_version() -> str:
    try:
        return f"primer3-py {primer3.__version__}"
    except Exception:
        return "primer3-py"



def write_csv(path: str, loci: list[Locus], tm: float) -> None:
    """One row per *ordered* pair -- primary and backup -- not one per locus.

    The earlier version emitted only `pairs[0]`, which meant the file everyone
    actually orders from did not contain half the oligos being ordered. `id`
    (`<locus>_p<rank>`) is the unique key, because `name` now repeats; the
    downstream checkers (verify_primers.py, blast_specificity.py) key on `id`.

    `backup_why` is repeated on both rows of a locus on purpose: a spreadsheet
    row gets copied out of its file, and a non-independent backup must not
    travel without the sentence saying so.
    """
    cols = ["id", "name", "gene", "rank", "pair", "role", "independent",
            "exon", "aa", "chrom", "cut", "strand",
            "tier", "tier_why", "fwd", "rev", "len_f", "len_r",
            "tm_f", "tm_r", "dtm", "gc_f", "gc_r", "product_bp",
            "cut_from_left", "cut_from_right", "snp_fwd", "snp_rev", "snp_3prime",
            "hetero_dg", "ispcr_products", "ispcr_off_target", "ispcr_hits",
            "penalty",
            "fwd_with_overhang", "rev_with_overhang", "len_f_oh", "len_r_oh",
            "backup_why", "notes"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for loc in loci:
            if not loc.pairs:
                w.writerow(dict(id=loc.name, name=loc.name, gene=loc.gene,
                                rank=loc.rank,
                                exon=loc.exon, aa=loc.aa, chrom=loc.chrom,
                                cut=loc.cut, strand=loc.strand, tier="FAILED",
                                tier_why="no pair at any tier",
                                notes="; ".join(loc.notes)))
                continue
            by_i = {p["i"]: p for p in loc.pairs}
            for n, i in enumerate(loc.recommended):
                p = by_i[i]
                w.writerow(dict(
                    id=f"{loc.name}_p{i}",
                    name=loc.name, gene=loc.gene, rank=loc.rank,
                    pair=i, role="primary" if n == 0 else "backup",
                    independent=("" if n == 0
                                 else ("yes" if loc.rec_independent else "no")),
                    exon=loc.exon,
                    aa=loc.aa, chrom=loc.chrom, cut=loc.cut, strand=loc.strand,
                    tier=loc.tier["name"] if loc.snp_ok
                         else "UNVERIFIED-" + loc.tier["name"],
                    tier_why=loc.tier["why"] if loc.snp_ok
                             else "variant lookup failed -- NOT SNP-screened",
                    fwd=p["fwd"], rev=p["rev"],
                    len_f=len(p["fwd"]), len_r=len(p["rev"]),
                    tm_f=p["tm_f"], tm_r=p["tm_r"], dtm=p["dtm"],
                    gc_f=p["gc_f"], gc_r=p["gc_r"], product_bp=p["size"],
                    cut_from_left=p["cut_from_left"],
                    cut_from_right=p["cut_from_right"],
                    snp_fwd=p["snp_fwd"], snp_rev=p["snp_rev"],
                    snp_3prime=p["snp_3prime"],
                    hetero_dg=p.get("hetero_dg"),
                    ispcr_products=p.get("ispcr_n", ""),
                    ispcr_off_target=("" if p.get("ispcr_off") is None
                                      else p["ispcr_off"]),
                    ispcr_hits=" | ".join(p.get("ispcr_hits", [])[:4]),
                    penalty=p["penalty"],
                    fwd_with_overhang=OVERHANG_F + p["fwd"],
                    rev_with_overhang=OVERHANG_R + p["rev"],
                    len_f_oh=len(OVERHANG_F) + len(p["fwd"]),
                    len_r_oh=len(OVERHANG_R) + len(p["rev"]),
                    backup_why=loc.rec_why,
                    notes="; ".join(loc.notes)))


def write_order_sheet(path: str, loci: list[Locus], overhang: bool) -> None:
    """A flat, paste-friendly table: gene, guide, which primer, sequence.

    Four columns and nothing else. Scale and purification are identical for
    every oligo here (25 nmol, desalted), so repeating them on every row was
    noise -- they belong on the supplier's form, not in a list someone has to
    read. Anyone receiving this can paste it straight into an order or a mail.

    Writes the *recommended* pairs -- primary and backup -- not just the
    primary, since that is what the report preselects and what gets ordered.
    A per-gene copy is written alongside the combined file so one gene can be
    handed over on its own.
    """
    cols = ["Gene", "Guide", "Primer", "Sequence"]
    rows: list[list[str]] = []
    for loc in loci:
        if not loc.pairs:
            continue
        by_i = {p["i"]: p for p in loc.pairs}
        for i in loc.recommended:
            p = by_i[i]
            f = (OVERHANG_F + p["fwd"]) if overhang else p["fwd"]
            r = (OVERHANG_R + p["rev"]) if overhang else p["rev"]
            rows.append([loc.gene, f"g{loc.rank}", f"pair {i} forward", f])
            rows.append([loc.gene, f"g{loc.rank}", f"pair {i} reverse", r])

    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        w.writerows(rows)

    base, ext = os.path.splitext(path)
    for gene in dict.fromkeys(r[0] for r in rows):
        gpath = f"{base}_{gene}{ext}"
        with open(gpath, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            w.writerows(r for r in rows if r[0] == gene)


# ========================================================================= main
def main() -> None:
    global SORT_MODE, SPEC_META, PRODUCT_OPT
    ap = argparse.ArgumentParser(
        description="Tm-harmonised, SNP-aware genotyping primers for CRISPR cut sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--from", dest="sources", action="append", required=True,
                    metavar="CSV:RANKS",
                    help="pipeline CSV and the ranks to design for, e.g. "
                         "out/yap1_danRer11.csv:5,16,17 (repeatable)")
    ap.add_argument("--build", default="danRer11")
    ap.add_argument("-o", "--out", default="primers_out")
    ap.add_argument("--max-product", type=int, default=233,
                    help="hard cap on gene-specific amplicon length")
    ap.add_argument("--min-product", type=int, default=175,
                    help="hard lower bound on product size (default 175). The "
                         "MiSeq protocol wants ~200 bp; below ~175 the cut "
                         "site crowds the read ends")
    ap.add_argument("--opt-product", type=int, default=200,
                    help="product size to aim for (default 200, the MiSeq "
                         "amplicon optimum; ~267 bp with the Illumina tails)")
    ap.add_argument("--min-flank", type=int, default=40,
                    help="minimum bp from the cut to either product end")
    ap.add_argument("--tm", type=float, default=None,
                    help="fix the shared Tm instead of scanning for it")
    ap.add_argument("--tm-lo", type=float, default=57.0)
    ap.add_argument("--tm-hi", type=float, default=65.0)
    ap.add_argument("--tm-step", type=float, default=0.25)
    ap.add_argument("--flank-window", type=int, default=450,
                    help="bp fetched either side of the cut")
    ap.add_argument("--no-snpfisher", action="store_true",
                    help="skip the Zv9/SNPfisher strain query (Ensembl only)")
    ap.add_argument("--ispcr-bin", default=None,
                    help="path to UCSC's isPcr binary; enables genome-wide "
                         "specificity checking (needs --genome)")
    ap.add_argument("--genome", default=None,
                    help="local danRer11 .2bit or .fa for isPcr")
    ap.add_argument("--overhang-in-order-sheet", action="store_true",
                    help="prepend the Illumina overhangs in order_sheet.csv")
    ap.add_argument("--num-return", type=int, default=40,
                    help="pairs to ask Primer3 for per locus (default 40). "
                         "Deeper than it looks like it needs to be: the search "
                         "for an *independent* backup is what needs the depth, "
                         "and measured on a twelve-guide set 40 gives 8/12 loci a "
                         "genuine independent backup where 20 gives 7/12. "
                         "60 adds nothing further")
    ap.add_argument("--max-candidates", type=int, default=8,
                    help="how many ranked pairs per locus to keep in "
                         "candidates.json for the report (default 8)")
    ap.add_argument("--primer-min", type=int, default=SIZE["min"],
                    help=f"minimum primer length (default {SIZE['min']})")
    ap.add_argument("--primer-opt", type=int, default=SIZE["opt"],
                    help="length Primer3 aims for. Pushed near the maximum "
                         f"(default {SIZE['opt']}) because a longer exact match "
                         "is rarer in the genome; see the length-policy note")
    ap.add_argument("--primer-max", type=int, default=SIZE["max"],
                    help=f"maximum primer length (default {SIZE['max']})")
    ap.add_argument("--specificity", default=None, metavar="JSON",
                    help="specificity.json from ispcr_scan.py. Pairs are then "
                         "ranked specific-first and the result is carried into "
                         "the CSV, the order sheet and the report")
    ap.add_argument("--sort", choices=["length", "penalty"], default=SORT_MODE,
                    help="shortlist ordering: 'length' prefers longer oligos "
                         "within a 0.1 C dTm bucket (default), 'penalty' is "
                         "Primer3's own ranking")
    ap.add_argument("--mv", type=float, default=50.0, help="monovalent cation mM")
    ap.add_argument("--dv", type=float, default=1.5, help="divalent cation (Mg2+) mM")
    ap.add_argument("--dntp", type=float, default=0.6, help="dNTP mM")
    ap.add_argument("--dna", type=float, default=50.0, help="primer conc nM")
    ap.add_argument("--delay", type=float, default=0.4)
    args = ap.parse_args()

    conc = dict(mv=args.mv, dv=args.dv, dntp=args.dntp, dna=args.dna)
    SIZE.update(min=args.primer_min, opt=args.primer_opt, max=args.primer_max)
    SORT_MODE = args.sort
    if not SIZE["min"] <= SIZE["opt"] <= SIZE["max"]:
        sys.exit(f"primer lengths must satisfy min <= opt <= max, got "
                 f"{SIZE['min']}/{SIZE['opt']}/{SIZE['max']}")
    PRODUCT_OPT = args.opt_product
    if not args.min_product <= args.opt_product <= args.max_product:
        sys.exit(f"products must satisfy min <= opt <= max, got "
                 f"{args.min_product}/{args.opt_product}/{args.max_product}")
    for t in TIERS:
        t["pmax"] = args.max_product
        t["pmin"] = max(60, args.min_product - t["pmin_drop"])
        # the optimum is clamped into each tier's own range, so a fallback tier
        # that allows shorter products still aims as high as it is allowed to
        t["popt"] = min(max(args.opt_product, t["pmin"]), t["pmax"])
        t["flank"] = max(15, args.min_flank - t["flank_drop"])

    loci: list[Locus] = []
    for s in args.sources:
        loci.extend(read_loci(s))
    os.makedirs(args.out, exist_ok=True)
    log(f"[1/5] {len(loci)} loci: " + ", ".join(l.name for l in loci))

    log("[2/5] fetching templates from UCSC ...")
    for loc in loci:
        fetch_template(loc, args.build, args.flank_window)
        time.sleep(args.delay)

    log("[3/5] collecting variants ...")
    for loc in loci:
        try:
            for v in ensembl_variants(loc):
                loc.snps.setdefault(v["off"], []).append(v)
        except RuntimeError as e:
            loc.snp_ok = False
            loc.notes.append(f"VARIANT LOOKUP FAILED -- primers are NOT "
                             f"SNP-screened for this locus ({e})")
            log(f"      !! {loc.name}: {e}")
        time.sleep(args.delay)
    if not args.no_snpfisher:
        snpfisher_variants(loci, args.build, os.path.join(
            os.path.dirname(os.path.abspath(__file__)), ".cache"))
    for loc in loci:
        log(f"      {loc.name}: {len(loc.snps)} variant positions in "
            f"{len(loc.template)} bp")

    if args.tm is None:
        log(f"[4/5] scanning Tm {args.tm_lo}-{args.tm_hi} C ...")
        scan = scan_tm(loci, args.tm_lo, args.tm_hi, args.tm_step, conc)
        tm = pick_tm(scan)
        best = max(r["n_ok"] for r in scan)
        log(f"      best shared Tm = {tm} C ({best}/{len(loci)} loci at strict tier)")
        for r in scan:
            if r["n_ok"] >= max(1, best - 2):
                log(f"        {r['tm']:>5} C  {r['n_ok']:>2}/{len(loci)}  "
                    f"penalty {r['mean_penalty']}")
        with open(os.path.join(args.out, "tm_scan.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["tm", "n_ok", "mean_penalty"])
            w.writeheader()
            w.writerows(scan)
    else:
        tm = args.tm
        scan = []
        log(f"[4/5] shared Tm fixed at {tm} C")

    # Specificity, if it has been measured. Keyed on the oligo pair rather
    # than on a rank, because ranks are exactly what this run is about to
    # change -- a pair keeps its identity across re-designs, its position does
    # not.
    spec_map: dict[tuple[str, str], dict] = {}
    spec_meta = None
    if args.specificity:
        with open(args.specificity, encoding="utf-8") as fh:
            sd = json.load(fh)
        spec_meta = {k: sd.get(k) for k in
                     ("genome", "min_perfect", "max_mismatch", "max_product",
                      "generated")}
        for r in sd.get("pairs", []):
            spec_map[(r["fwd"], r["rev"])] = r
        SPEC_META = spec_meta
        log(f"      specificity loaded for {len(spec_map)} pairs "
            f"({spec_meta.get('genome')}, 3' seed {spec_meta.get('min_perfect')} nt, "
            f"<= {spec_meta.get('max_mismatch')} mismatches)")

    log("[5/5] designing, walking the relaxation ladder where needed ...")
    for loc in loci:
        for tier in TIERS:
            pairs = design_locus(loc, tm, tier, conc, args.num_return)
            if pairs:
                for pr in pairs:
                    sp = spec_map.get((pr["fwd"], pr["rev"]))
                    if sp:
                        pr["ispcr_n"] = sp.get("products")
                        pr["ispcr_off"] = sp.get("off_target")
                        pr["ispcr_verdict"] = sp.get("verdict")
                        pr["ispcr_hits"] = [f"{a['chrom']}:{a['start']}-{a['end']}"
                                            f"({a['size']}bp,{a['mm']}mm)"
                                            for a in (sp.get("off") or [])[:4]]
                pairs.sort(key=pair_sort_key)
                for n, pr in enumerate(pairs, 1):
                    pr["i"] = n          # rank, used by the order sheet and report
                loc.pairs, loc.tier = pairs, tier
                break
        flag = loc.tier["name"] if loc.tier else "FAILED"
        if loc.pairs:
            p = loc.pairs[0]
            log(f"      {loc.name:<12} {flag:<11} {p['size']:>3} bp  "
                f"Tm {p['tm_f']}/{p['tm_r']} (d{p['dtm']})  "
                f"cut {p['cut_from_left']}/{p['cut_from_right']}  "
                f"SNPs {p['snp_fwd']}+{p['snp_rev']}")
        else:
            log(f"      {loc.name:<12} FAILED -- no pair at any tier")

    for loc in loci:
        (loc.recommended, loc.rec_why,
         loc.rec_independent) = recommend_pairs(loc.pairs)

    if args.ispcr_bin and args.genome:
        log("      in-silico PCR (local isPcr) ...")
        req = [(loc.name, loc.pairs[0]["fwd"], loc.pairs[0]["rev"])
               for loc in loci if loc.pairs]
        hits = ispcr_local(args.ispcr_bin, args.genome, req)
        for loc in loci:
            if not loc.pairs:
                continue
            h = hits.get(loc.name)
            if h is None:
                continue
            p = loc.pairs[0]
            p["ispcr_n"], p["ispcr_hits"] = len(h), h
            if len(h) == 0:
                # Cannot be a real answer: the pair came off this assembly.
                loc.notes.append("isPcr found 0 products, impossible for a pair "
                                 "designed on this assembly -- check the genome "
                                 "file; specificity NOT verified")
            elif len(h) > 1:
                loc.notes.append(f"isPcr: {len(h)} products -- check specificity")
            log(f"      {loc.name:<12} {len(h)} product(s)")
    elif spec_map:
        n_known = n_clean = 0
        for loc in loci:
            for pr in loc.pairs:
                if pr.get("ispcr_off") is None:
                    continue
                n_known += 1
                n_clean += pr["ispcr_off"] == 0
            rec_off = [pr.get("ispcr_off") for pr in loc.pairs
                       if pr["i"] in loc.recommended]
            if any(o is None for o in rec_off):
                loc.notes.append("specificity not checked for every ordered pair")
            elif any(o for o in rec_off):
                loc.notes.append(f"in-silico PCR: {max(rec_off)} off-target "
                                 f"product(s) on an ordered pair -- see "
                                 f"ispcr_hits")
            else:
                loc.notes.append(f"in-silico PCR on {spec_meta.get('genome')}: "
                                 f"single product")
        log(f"      specificity: in-silico PCR on {spec_meta.get('genome')} -- "
            f"{n_clean}/{n_known} candidate pairs predict exactly one product")
        # Overwrite, never leave behind, the "paste these into Primer-BLAST"
        # file from an earlier run: a stale to-do list that says specificity is
        # unchecked is worse than no file, because it is read as current.
        with open(os.path.join(args.out, "specificity_todo.txt"), "w",
                  encoding="utf-8") as fh:
            fh.write("Specificity WAS checked for this design -- nothing to do "
                     "by hand.\n\n"
                     f"In-silico PCR against {spec_meta.get('genome')} "
                     f"(ispcr_scan.py, {spec_meta.get('generated')}): a site "
                     f"must match the primer's 3' "
                     f"{spec_meta.get('min_perfect')} nt exactly, with at most "
                     f"{spec_meta.get('max_mismatch')} mismatches 5' of that, "
                     f"and products up to {spec_meta.get('max_product')} bp are "
                     f"paired up.\n\n"
                     "Per-pair results: primers_out/specificity.json, the "
                     "ispcr_* columns of primers.csv, and the In-silico PCR "
                     "column of the report.\n\n")
            flagged = [(loc, pr) for loc in loci for pr in loc.pairs
                       if pr["i"] in loc.recommended and pr.get("ispcr_off")]
            if flagged:
                fh.write("Ordered pairs with a predicted off-target product:\n")
                for loc, pr in flagged:
                    fh.write(f"  {loc.name}_p{pr['i']}  "
                             f"{pr['ispcr_off']} off-target: "
                             f"{'; '.join(pr.get('ispcr_hits') or [])}\n")
            else:
                fh.write("Every ordered pair predicts a single product.\n")
            fh.write("\nThe reference assembly is not your fish: a test PCR is "
                     "still the real control.\n")
    else:
        log("      specificity: NOT CHECKED (no --ispcr-bin/--genome). "
            "UCSC's web in-silico PCR is behind a bot check and must not be "
            "scraped; run Primer-BLAST by hand or supply a local genome.")
        for loc in loci:
            loc.notes.append("specificity not checked")
        with open(os.path.join(args.out, "specificity_todo.txt"), "w") as fh:  # noqa: E501
            fh.write("Paste each pair into NCBI Primer-BLAST "
                     "(https://www.ncbi.nlm.nih.gov/tools/primer-blast/),\n"
                     "organism Danio rerio, and confirm a single product.\n\n")
            for loc in loci:
                by_i = {p["i"]: p for p in loc.pairs}
                for i in loc.recommended:
                    p = by_i[i]
                    fh.write(f"{loc.name}_p{i}\n  F  {p['fwd']}\n"
                             f"  R  {p['rev']}\n\n")

    csv_path = os.path.join(args.out, "primers.csv")
    write_csv(csv_path, loci, tm)
    write_order_sheet(os.path.join(args.out, "order_sheet.csv"), loci,
                      args.overhang_in_order_sheet)
    cand_path = os.path.join(args.out, "candidates.json")
    write_candidates(cand_path, loci, tm, args.max_candidates, args, scan)

    ok = sum(1 for l in loci if l.pairs)
    strict = sum(1 for l in loci if l.tier and l.tier["name"] == "strict")
    log(f"\ndone: {ok}/{len(loci)} loci have primers, {strict} at the strict tier, "
        f"shared Tm {tm} C")
    log(f"  {csv_path}")
    log(f"  {os.path.join(args.out, 'order_sheet.csv')}")


if __name__ == "__main__":
    main()
