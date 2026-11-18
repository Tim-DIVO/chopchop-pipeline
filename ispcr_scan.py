#!/usr/bin/env python3
"""
ispcr_scan.py -- genome-wide in-silico PCR, locally, with no web service.

Why
---
Every specificity route this project tried is blocked or unreliable:

* UCSC's web In-Silico PCR (`hgPcr`) sits behind a Cloudflare bot check that
  points automated clients at the downloads FAQ. Getting past it means lying
  about the User-Agent, so the tool does not.
* UCSC's `isPcr` binary is the sanctioned answer, but it is distributed as a
  Linux/macOS executable and this workstation is Windows.
* NCBI's BLAST URL API works, but the public queue is the bottleneck: a
  2-primer search submitted with an RTOE estimate of 5 s sat `WAITING` for over
  20 minutes (measured 2026-11-16), and a 24-pair search is not faster for
  being bigger. It is a fine background check, not something to design against.

So this does the same job directly against a local `.2bit`: find where each
primer can actually prime, pair those sites up, and report the products.

The algorithm, and why it is this one
-------------------------------------
It is isPcr's, deliberately, because isPcr's model of priming is the right one:

1. **Perfect 3' seed.** A site must match the primer's last `--min-perfect`
   bases (default 15) exactly. A mismatch inside the last few bases of a primer
   stops extension; a mismatch at the 5' end usually does not. Matching on the
   3' end is what separates "a sequence that looks a bit like the primer" from
   "a site that will amplify".
2. **Tolerant 5' tail.** 5' of the seed, up to `--max-mismatch` mismatches
   (default 3) are allowed. No gaps: a bulged primer that still primes is rare
   enough that modelling it would add far more false positives than it removes.
3. **Pair up.** A product is a plus-strand site and a minus-strand site on the
   same sequence, facing each other, within `--max-product` bp (default 3000,
   well past anything a 30 s extension would give).

Counting raw BLAST-style hits without step 1 over-reports off-targets wildly,
which is the trap this avoids.

What it is not
--------------
Not a wet-lab control. It searches the reference assembly: a product unique in
danRer11 can still have a partner in a line carrying an insertion, and nothing
here predicts how a PCR behaves at the bench.

Usage
-----
  python ispcr_scan.py --genome .cache/genome/danRer11.2bit \\
      --candidates primers_out/candidates.json \\
      -o primers_out/specificity.json

  # or straight from primers.csv (ordered pairs only)
  python ispcr_scan.py --genome ... --primers primers_out/primers.csv -o ...
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import struct
import sys
import time
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("numpy is required:  pip install numpy")

COMP = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")


def seq_class(name: str) -> str:
    """primary | unplaced | alt -- and the distinction decides correctness.

    danRer11 is 26 primary sequences, 967 unplaced scaffolds, and **930 alt
    contigs**. An `_alt` (or `_fix`) contig is an *alternative representation of
    a region that is already in the primary assembly* -- the same piece of
    genome, written down twice. A fish carries one version of it, so a product
    predicted on an alt contig is not a second band: it is the intended product
    seen through the assembly's second copy.

    Counting them as off-targets is not a harmless excess of caution. On this
    very run it flagged **all twelve** ordered yap1 pairs -- every candidate at
    every yap1 locus, because `chr18_KZ115577v1_alt` covers the whole gene -- so
    ranking on it would have demoted good primers for a reason that does not
    exist at the bench, and picked worse ones in their place.

    Unplaced scaffolds (`chrUn_*`, `*_random`) are the opposite case: real
    sequence that simply has no coordinate yet, so a product there is a real
    second product and counts.
    """
    if re.fullmatch(r"chr[0-9]+|chrM", name):
        return "primary"
    if name.endswith("_alt") or name.endswith("_fix"):
        return "alt"
    return "unplaced"


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


# ============================================================ .2bit reading
class TwoBit:
    """Minimal reader for UCSC's .2bit format.

    Format (little- or big-endian, signature 0x1A412743): a header, then an
    index of name -> file offset, then per sequence: dnaSize, a list of N
    blocks, a list of mask blocks, then DNA packed 4 bases per byte in the
    order T,C,A,G (00,01,10,11).

    Mask blocks (soft-masking = repeats) are read but NOT applied: a primer
    sitting in a repeat is exactly the case worth reporting, so masking it out
    would hide the problem instead of finding it.
    """

    SIG = 0x1A412743

    def __init__(self, path: str):
        self.fh = open(path, "rb")
        sig = struct.unpack("<I", self.fh.read(4))[0]
        if sig == self.SIG:
            self.end = "<"
        else:
            self.fh.seek(0)
            sig = struct.unpack(">I", self.fh.read(4))[0]
            if sig != self.SIG:
                raise ValueError(f"{path} is not a .2bit file")
            self.end = ">"
        ver, count, _ = struct.unpack(self.end + "III", self.fh.read(12))
        self.index: dict[str, int] = {}
        for _ in range(count):
            nlen = self.fh.read(1)[0]
            name = self.fh.read(nlen).decode()
            off = struct.unpack(self.end + "I", self.fh.read(4))[0]
            self.index[name] = off
        # byte -> 4 bases, the 2bit code order
        table = np.frombuffer(b"TCAG", dtype="S1")
        codes = np.arange(256, dtype=np.uint8)
        self._lut = np.stack([table[(codes >> s) & 3] for s in (6, 4, 2, 0)], axis=1)

    def names(self) -> list[str]:
        return list(self.index)

    def sequence(self, name: str) -> bytes:
        """Uppercase ACGT bytes, with N blocks written as N."""
        self.fh.seek(self.index[name])
        (dna_size,) = struct.unpack(self.end + "I", self.fh.read(4))
        (n_count,) = struct.unpack(self.end + "I", self.fh.read(4))
        n_starts = np.frombuffer(self.fh.read(4 * n_count), dtype=self.end + "u4")
        n_sizes = np.frombuffer(self.fh.read(4 * n_count), dtype=self.end + "u4")
        (m_count,) = struct.unpack(self.end + "I", self.fh.read(4))
        self.fh.seek(8 * m_count, os.SEEK_CUR)      # mask blocks: skipped, see above
        self.fh.seek(4, os.SEEK_CUR)                # reserved
        packed = np.frombuffer(self.fh.read((dna_size + 3) // 4), dtype=np.uint8)
        seq = self._lut[packed].reshape(-1)[:dna_size]
        for s, ln in zip(n_starts.tolist(), n_sizes.tolist()):
            seq[s:s + ln] = b"N"
        return seq.tobytes()


# ============================================================ site finding
def find_sites(seq: bytes, primers: list[tuple[str, str]], min_perfect: int,
               max_mm: int) -> dict[str, list[tuple[int, str, int]]]:
    """All primable sites for each primer in one sequence.

    Returns name -> [(genomic 0-based start of the primer's 5' end, strand,
    mismatches)]. Plus-strand sites are found by searching the seed itself;
    minus-strand sites by searching the seed's reverse complement, which is
    what the primer would anneal to when priming leftwards.
    """
    out: dict[str, list[tuple[int, str, int]]] = defaultdict(list)
    for name, p in primers:
        pb = p.encode().upper()
        L = len(pb)
        k = min(min_perfect, L)
        seed = pb[L - k:]                       # the 3' end
        tail = pb[:L - k]                       # 5' of the seed

        # --- plus strand: ...tail|seed-> , primer extends rightwards
        i = seq.find(seed)
        while i >= 0:
            start = i - len(tail)
            if start >= 0:
                mm = _mismatches(seq, start, tail)
                if mm <= max_mm:
                    out[name].append((start, "+", mm))
            i = seq.find(seed, i + 1)

        # --- minus strand: <-seed|tail... on the plus strand, reverse complemented
        rseed = seed.translate(COMP)[::-1]
        rtail = tail.translate(COMP)[::-1]
        i = seq.find(rseed)
        while i >= 0:
            tstart = i + len(rseed)
            if tstart + len(rtail) <= len(seq):
                mm = _mismatches(seq, tstart, rtail)
                if mm <= max_mm:
                    # report the primer's 3'-most genomic coordinate first:
                    # the site spans [i, i+len(primer))
                    out[name].append((i, "-", mm))
            i = seq.find(rseed, i + 1)
    return out


def _mismatches(seq: bytes, start: int, sub: bytes) -> int:
    if not sub:
        return 0
    window = seq[start:start + len(sub)]
    if len(window) < len(sub):
        return len(sub)          # runs off the end: unusable
    return sum(1 for a, b in zip(window, sub) if a != b)


def products(f_sites, r_sites, flen: int, rlen: int, max_product: int):
    """Pair plus- and minus-strand sites into predicted amplicons.

    Either primer can play either role at an off-target site, so both
    assignments are tried; identical predictions are collapsed.
    """
    out = []
    for a_sites, a_len, b_sites, b_len in ((f_sites, flen, r_sites, rlen),
                                           (r_sites, rlen, f_sites, flen)):
        plus = [s for s in a_sites if s[1] == "+"]
        minus = [s for s in b_sites if s[1] == "-"]
        for ps, _, pmm in plus:
            for ms, _, mmm in minus:
                end = ms + b_len               # one past the product's last base
                size = end - ps
                # The two primer sites must not overlap and must face each
                # other: a minus site that starts before the plus site ends is
                # not a product, it is two primers pointing the wrong way, and
                # `size > 0` alone lets a few of those through.
                if ms >= ps + a_len and a_len + b_len <= size <= max_product:
                    out.append((ps, end, size, pmm + mmm))
    seen, uniq = set(), []
    for p in sorted(out, key=lambda x: (x[3], x[2])):
        if p[:2] in seen:
            continue
        seen.add(p[:2])
        uniq.append(p)
    return uniq


# ============================================================ driver
def load_pairs(args) -> list[dict]:
    """Every (locus, pair) to check, from candidates.json or primers.csv."""
    pairs = []
    if args.candidates:
        d = json.load(open(args.candidates, encoding="utf-8"))
        for L in d["loci"]:
            for p in L["pairs"]:
                pairs.append(dict(locus=L["name"], gene=L["gene"], i=p["i"],
                                  fwd=p["fwd"], rev=p["rev"], size=p["size"],
                                  chrom=L["chrom"],
                                  recommended=p["i"] in (L.get("recommended") or [])))
    else:
        for r in csv.DictReader(open(args.primers, newline="")):
            if not r.get("fwd"):
                continue
            pairs.append(dict(locus=r["name"], gene=r.get("gene", ""),
                              i=int(r.get("pair") or 1), fwd=r["fwd"], rev=r["rev"],
                              size=int(r["product_bp"]), chrom=r.get("chrom", ""),
                              recommended=True))
    return pairs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help="UCSC .2bit")
    ap.add_argument("--candidates", default=None, help="candidates.json (all pairs)")
    ap.add_argument("--primers", default=None, help="primers.csv (ordered pairs only)")
    ap.add_argument("-o", "--out", default=None, help="write results JSON here")
    ap.add_argument("--min-perfect", type=int, default=15,
                    help="bases at the 3' end that must match exactly (isPcr's "
                         "default, and the parameter that matters most)")
    ap.add_argument("--max-mismatch", type=int, default=3,
                    help="mismatches tolerated 5' of the seed (default 3)")
    ap.add_argument("--max-product", type=int, default=3000)
    ap.add_argument("--only-main", action="store_true",
                    help="skip alt/random/unplaced scaffolds entirely")
    ap.add_argument("--count-alt", action="store_true",
                    help="count products on alt/fix contigs as off-targets. "
                         "Off by default: an alt contig is the same locus "
                         "written twice, not a second site (see seq_class)")
    args = ap.parse_args()
    if not (args.candidates or args.primers):
        sys.exit("need --candidates or --primers")

    pairs = load_pairs(args)
    primers: dict[str, str] = {}
    for p in pairs:
        primers[f"{p['locus']}_p{p['i']}_F"] = p["fwd"]
        primers[f"{p['locus']}_p{p['i']}_R"] = p["rev"]
    # the same oligo often appears in several pairs; search each sequence once
    uniq_seq: dict[str, str] = {}
    for n, s in primers.items():
        uniq_seq.setdefault(s, n)
    plist = [(n, s) for s, n in uniq_seq.items()]
    print(f"{len(pairs)} pairs, {len(primers)} primer slots, "
          f"{len(plist)} distinct oligos")

    tb = TwoBit(args.genome)
    names = tb.names()
    if args.only_main:
        names = [n for n in names if "_" not in n]
    print(f"genome: {os.path.basename(args.genome)}, {len(names)} sequences")

    sites: dict[str, list] = defaultdict(list)      # oligo seq -> [(chrom,pos,strand,mm)]
    t0 = time.time()
    for n, name in enumerate(names, 1):
        seq = tb.sequence(name)
        found = find_sites(seq, plist, args.min_perfect, args.max_mismatch)
        for qn, hits in found.items():
            s = primers[qn]
            for pos, strand, mm in hits:
                sites[s].append((name, pos, strand, mm))
        if n % 25 == 0 or n == len(names):
            print(f"  {n}/{len(names)} sequences, {time.time() - t0:.0f}s, "
                  f"{sum(len(v) for v in sites.values())} sites", flush=True)
        del seq

    # ---- pair the sites up, per primer pair
    results = []
    for p in pairs:
        f, r = sites.get(p["fwd"], []), sites.get(p["rev"], [])
        by_chrom_f, by_chrom_r = defaultdict(list), defaultdict(list)
        for c, pos, strand, mm in f:
            by_chrom_f[c].append((pos, strand, mm))
        for c, pos, strand, mm in r:
            by_chrom_r[c].append((pos, strand, mm))
        amps = []
        for c in set(by_chrom_f) | set(by_chrom_r):
            for lo, hi, size, mm in products(by_chrom_f.get(c, []),
                                             by_chrom_r.get(c, []),
                                             len(p["fwd"]), len(p["rev"]),
                                             args.max_product):
                amps.append(dict(chrom=c, start=lo + 1, end=hi, size=size, mm=mm))
        amps.sort(key=lambda a: (a["mm"], a["size"]))
        for a in amps:
            a["kind"] = seq_class(a["chrom"])
        intended = [a for a in amps if a["size"] == p["size"]
                    and (not p["chrom"] or a["chrom"] == p["chrom"])]
        rest = [a for a in amps if a not in intended]
        # alt-haplotype copies are reported, never counted (unless asked)
        alt = [] if args.count_alt else [a for a in rest if a["kind"] == "alt"]
        extra = [a for a in rest if a not in alt]
        verdict = ("INTENDED PRODUCT NOT FOUND" if not intended else
                   "specific" if not extra else
                   f"{len(extra)} off-target product(s)")
        if verdict == "specific" and alt:
            verdict = f"specific ({len(alt)} alt-haplotype copy)"
        results.append(dict(locus=p["locus"], gene=p["gene"], pair=p["i"],
                            fwd=p["fwd"], rev=p["rev"], size=p["size"],
                            recommended=p["recommended"],
                            f_sites=len(f), r_sites=len(r),
                            products=len(amps), off_target=len(extra),
                            alt_copies=len(alt), verdict=verdict,
                            off=extra[:5], alt=alt[:3]))

    bad = [r for r in results if r["off_target"] or not r["verdict"].startswith("specific")]
    n_alt = sum(1 for r in results if r.get("alt_copies"))
    print(f"\n{len(results) - len(bad)}/{len(results)} candidate pairs are "
          f"predicted specific in the primary assembly.")
    if n_alt:
        print(f"{n_alt} of them also hit an alt-haplotype contig -- reported, "
              f"not counted (the same locus written twice; --count-alt to "
              f"include).")
    rec = [r for r in results if r["recommended"]]
    rbad = [r for r in rec if r["off_target"] or not r["verdict"].startswith("specific")]
    if rec:
        print(f"of the {len(rec)} currently recommended pairs: "
              f"{len(rec) - len(rbad)} specific, {len(rbad)} not")
    for r in rbad:
        print(f"  !! {r['locus']}_p{r['pair']}: {r['verdict']}")
        for a in r["off"][:3]:
            print(f"       {a['chrom']}:{a['start']}-{a['end']} "
                  f"({a['size']} bp, {a['mm']} mismatch)")

    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            json.dump({"genome": os.path.basename(args.genome),
                       "min_perfect": args.min_perfect,
                       "max_mismatch": args.max_mismatch,
                       "max_product": args.max_product,
                       "count_alt": args.count_alt,
                       "generated": time.strftime("%Y-%m-%d %H:%M"),
                       "pairs": results}, fh)
        print(f"\nwrote {args.out}")
    print("\nThe reference assembly is not your fish. A test PCR is still the "
          "real control.")


if __name__ == "__main__":
    main()
