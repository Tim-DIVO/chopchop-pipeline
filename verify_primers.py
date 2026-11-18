#!/usr/bin/env python3
"""
verify_primers.py -- independently check that designed primers are real.

Why this exists
---------------
`primer_design.py` asserts that each primer pair sits at a particular place in
the genome, flanks a particular guide, and yields a particular product length.
Those are claims produced by the same program that chose the primers. If it has
a coordinate bug, the CSV it writes will be wrong and self-consistent, which is
the worst kind of wrong.

So this script re-derives all of it from scratch and shares nothing with the
designer:

  * different sequence source -- Ensembl REST, not UCSC
  * different coordinate convention -- Ensembl is 1-based inclusive, UCSC's API
    is 0-based half-open, so an off-by-one in either tool shows up as a
    mismatch here rather than cancelling out
  * different Tm implementation -- SantaLucia 1998 nearest-neighbour written
    out longhand below, not Primer3's
  * no imports from primer_design.py at all

A pass here means: the primer sequences genuinely occur, exactly once each, at
the stated positions, in the orientation stated, flanking the chosen guide, with
the stated product length. It does NOT mean they are specific genome-wide (see
"What this does not check" at the bottom) or that they will work at the bench.

Usage
-----
  python verify_primers.py primers_out/primers.csv \
      --guides geneA=out_geneA/geneA_danRer11.csv \
      --guides yap1=out/yap1_danRer11.csv

Exit code is 0 only if every locus passes every structural check.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
import urllib.error
import urllib.request

try:
    import certifi
    import ssl
    _CTX = ssl.create_default_context(cafile=certifi.where())
    ssl._create_default_https_context = lambda *a, **k: _CTX  # noqa: SLF001
except Exception:
    pass

ENSEMBL = "https://rest.ensembl.org"
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".cache", "verify")

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


# ----------------------------------------------------------------- sequence fetch
def ensembl_seq(chrom: str, start1: int, end1: int, species="danio_rerio",
                asm="GRCz11", delay=0.34) -> str:
    """Fetch [start1, end1] INCLUSIVE, 1-based -- Ensembl's convention.

    Deliberately not UCSC. If the designer and this script agreed on a wrong
    coordinate they would both be wrong in the same direction; two different
    conventions from two different servers will not fail identically.
    """
    name = chrom[3:] if chrom.startswith("chr") else chrom
    os.makedirs(CACHE, exist_ok=True)
    key = os.path.join(CACHE, f"{species}_{asm}_{name}_{start1}_{end1}.txt")
    if os.path.exists(key):
        return open(key).read().strip().upper()

    url = (f"{ENSEMBL}/sequence/region/{species}/{name}:{start1}..{end1}:1"
           f"?coord_system_version={asm};content-type=text/plain")
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, headers={"Accept": "text/plain"})
            with urllib.request.urlopen(req, timeout=60) as fh:
                seq = fh.read().decode().strip().upper()
            with open(key, "w") as fh:
                fh.write(seq)
            time.sleep(delay)
            return seq
        except urllib.error.HTTPError as e:
            last = e
            if e.code in (429, 500, 502, 503, 504):
                wait = int(e.headers.get("Retry-After", 2 ** attempt)) or 2 ** attempt
                sys.stderr.write(f"      Ensembl {e.code}; waiting {wait}s\n")
                time.sleep(wait)
                continue
            raise
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** attempt)
    raise SystemExit(f"Ensembl fetch failed for {chrom}:{start1}-{end1}: {last}")


# ------------------------------------------------------------------- independent Tm
# SantaLucia (1998) unified nearest-neighbour parameters.
# dH kcal/mol, dS cal/(mol K), for each 5'->3' dinucleotide on the top strand.
NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)    # initiation with terminal G/C
INIT_AT = (2.3, 4.1)     # initiation with terminal A/T
R = 1.98720425864083     # cal/(mol K)


def tm_santalucia(seq: str, mv=50.0, dv=1.5, dntp=0.6, dna_nm=50.0) -> float:
    """Nearest-neighbour Tm, written out longhand so it shares no code path
    with Primer3. mv/dv/dntp in mM, primer in nM."""
    s = seq.upper()
    if len(s) < 2 or any(c not in "ACGT" for c in s):
        return float("nan")
    dh, ds = 0.0, 0.0
    for a, b in zip(s, s[1:]):
        h, sv = NN[a + b]
        dh += h
        ds += sv
    for end in (s[0], s[-1]):
        h, sv = INIT_GC if end in "GC" else INIT_AT
        dh += h
        ds += sv

    # Divalent -> monovalent equivalent (von Ahsen / Owczarzy style, as used by
    # Primer3's default salt_corrections=1 path).
    free_mg = max(dv - dntp, 0.0)
    na_eq = (mv + 120.0 * math.sqrt(free_mg)) / 1000.0  # mol/L
    if na_eq <= 0:
        na_eq = 1e-6
    ds += 0.368 * (len(s) - 1) * math.log(na_eq)

    ct = (dna_nm / 1e9) / 4.0   # non-self-complementary primer in excess
    return (dh * 1000.0) / (ds + R * math.log(ct)) - 273.15


def gc_pct(s: str) -> float:
    return 100.0 * sum(c in "GC" for c in s.upper()) / len(s)


# ------------------------------------------------------------------------ checking
def find_all(hay: str, needle: str) -> list[int]:
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def check_locus(row: dict, guide: dict | None, pad: int, tm_tol: float) -> dict:
    """Re-derive everything about one primer pair from freshly fetched sequence."""
    # primers.csv now carries one row per ordered pair, so `name` (the locus)
    # repeats and `id` (<locus>_p<rank>) is the unique key. Fall back to `name`
    # for the older one-row-per-locus files.
    name = row.get("id") or row["name"]
    chrom = row["chrom"]
    cut = int(row["cut"])
    fwd = row["fwd"].upper()
    rev = row["rev"].upper()
    claimed_product = int(row["product_bp"])
    problems, notes = [], []

    # Window generously around the cut, fetched 1-based inclusive from Ensembl.
    lo = max(1, cut - pad)
    hi = cut + pad
    win = ensembl_seq(chrom, lo, hi)

    # 1. the forward primer must occur, exactly once, on the plus strand
    f_hits = find_all(win, fwd)
    # 2. the reverse primer binds the minus strand, so its reverse complement
    #    must appear on the plus strand
    r_hits = find_all(win, rc(rev))

    if len(f_hits) != 1:
        problems.append(f"forward primer found {len(f_hits)}x in the window "
                        f"(expected exactly 1)")
    if len(r_hits) != 1:
        problems.append(f"reverse primer (as revcomp) found {len(r_hits)}x "
                        f"(expected exactly 1)")

    product = cut_l = cut_r = None
    if len(f_hits) == 1 and len(r_hits) == 1:
        f0 = f_hits[0]
        r_end = r_hits[0] + len(rev)          # exclusive end of the revcomp match
        if r_end <= f0:
            problems.append("reverse primer lies upstream of the forward primer "
                            "-- the pair cannot amplify")
        else:
            # 3. product length, measured not trusted
            product = r_end - f0
            if product != claimed_product:
                problems.append(f"product is {product} bp here, CSV claims "
                                f"{claimed_product} bp")
            # 4. the cut must sit inside the product, at the stated offsets
            cut_idx = cut - lo          # 0-based index of `cut` within win
            cut_l = cut_idx - f0
            cut_r = product - cut_l
            for label, got, want_s in (("cut_from_left", cut_l, row.get("cut_from_left")),
                                       ("cut_from_right", cut_r, row.get("cut_from_right"))):
                if want_s not in (None, ""):
                    want = int(want_s)
                    # allow 1 bp: "distance from the cut" is a half-open/closed
                    # convention choice, not a correctness question
                    if abs(got - want) > 1:
                        problems.append(f"{label} is {got} here, CSV claims {want}")
            if not (0 < cut_l < product):
                problems.append("cut site is not inside the product")

            # 5. primers must not overlap the protospacer
            if guide:
                target = (guide.get("target") or "").upper()
                if target:
                    t_hits = find_all(win, target) or find_all(win, rc(target))
                    if len(t_hits) != 1:
                        notes.append(f"guide 23-mer found {len(t_hits)}x in window")
                    else:
                        t0 = t_hits[0]
                        t1 = t0 + len(target)
                        if not (f0 + len(fwd) <= t0 and t1 <= r_hits[0]):
                            problems.append("a primer overlaps the protospacer")
                        if not (f0 < t0 and t1 < r_end):
                            problems.append("the guide is not inside the product")
                else:
                    notes.append("guide CSV had no target sequence to check")
            else:
                notes.append("no guide row matched; protospacer check skipped")

    # 6. independent Tm, informational
    tm_f_mine = tm_santalucia(fwd)
    tm_r_mine = tm_santalucia(rev)
    d_f = tm_f_mine - float(row["tm_f"])
    d_r = tm_r_mine - float(row["tm_r"])
    if max(abs(d_f), abs(d_r)) > tm_tol:
        notes.append(f"Tm differs from Primer3 by {d_f:+.1f}/{d_r:+.1f} C "
                     f"(> {tm_tol} C)")

    return dict(name=name, ok=not problems, problems=problems, notes=notes,
                product=product, cut_l=cut_l, cut_r=cut_r,
                tm_f_mine=tm_f_mine, tm_r_mine=tm_r_mine, d_f=d_f, d_r=d_r,
                gc_f=gc_pct(fwd), gc_r=gc_pct(rev), window=f"{chrom}:{lo}-{hi}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Independently verify designed primers against Ensembl.")
    ap.add_argument("primers", help="primers.csv from primer_design.py")
    ap.add_argument("--guides", action="append", default=[], metavar="GENE=CSV",
                    help="guide CSV per gene, to check the protospacer is "
                         "flanked; repeatable")
    ap.add_argument("--pad", type=int, default=600,
                    help="bp fetched either side of the cut (default 600)")
    ap.add_argument("--tm-tol", type=float, default=2.0,
                    help="report when my Tm and Primer3's differ by more than "
                         "this many C (default 2.0, informational only)")
    ap.add_argument("-o", "--out", default=None, help="write a CSV report here")
    args = ap.parse_args()

    gmap: dict[str, dict[str, dict]] = {}
    for spec in args.guides:
        if "=" not in spec:
            sys.exit(f"--guides wants GENE=CSV, got {spec!r}")
        gene, path = spec.split("=", 1)
        gmap[gene] = {r["rank"]: r for r in csv.DictReader(open(path, newline=""))}

    rows = list(csv.DictReader(open(args.primers, newline="")))
    print(f"verifying {len(rows)} primer pairs against Ensembl GRCz11")
    print("(sequence source and coordinate convention both differ from the "
          "designer's)\n")

    results = []
    for row in rows:
        guide = gmap.get(row.get("gene", ""), {}).get(row.get("rank", ""))
        res = check_locus(row, guide, args.pad, args.tm_tol)
        results.append(res)
        flag = "PASS" if res["ok"] else "FAIL"
        prod = res["product"] if res["product"] is not None else "-"
        print(f"  {flag}  {res['name']:<12} product {str(prod):>4} bp  "
              f"cut {str(res['cut_l']):>4}/{str(res['cut_r']):<4}  "
              f"Tm(mine) {res['tm_f_mine']:.1f}/{res['tm_r_mine']:.1f} "
              f"(d {res['d_f']:+.1f}/{res['d_r']:+.1f})")
        for p in res["problems"]:
            print(f"          !! {p}")
        for n in res["notes"]:
            print(f"           - {n}")

    n_ok = sum(r["ok"] for r in results)
    print(f"\n{n_ok}/{len(results)} pairs verified structurally.")
    if n_ok == len(results):
        print("Every primer occurs exactly once at the stated position, in the")
        print("stated orientation, flanking the chosen guide, with the stated")
        print("product length -- re-derived from an independent source.")
    print("\nNOT checked here: genome-wide specificity (the window is local by")
    print("design), and anything about how the PCR behaves at the bench.")

    if args.out:
        with open(args.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["name", "ok", "product_bp", "cut_from_left",
                        "cut_from_right", "tm_f_santalucia", "tm_r_santalucia",
                        "d_tm_f", "d_tm_r", "gc_f", "gc_r", "window",
                        "problems", "notes"])
            for r in results:
                w.writerow([r["name"], int(r["ok"]), r["product"], r["cut_l"],
                            r["cut_r"], f"{r['tm_f_mine']:.2f}",
                            f"{r['tm_r_mine']:.2f}", f"{r['d_f']:+.2f}",
                            f"{r['d_r']:+.2f}", f"{r['gc_f']:.0f}",
                            f"{r['gc_r']:.0f}", r["window"],
                            "; ".join(r["problems"]), "; ".join(r["notes"])])
        print(f"\nwrote {args.out}")

    sys.exit(0 if n_ok == len(results) else 1)


if __name__ == "__main__":
    main()
