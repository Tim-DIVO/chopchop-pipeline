#!/usr/bin/env python3
"""
inject_primers.py -- re-point an existing report at the current candidates.json.

Why this exists
---------------
`chopchop_pipeline.py` injects the primer shortlist when it writes a report, by
substituting the `/*__PRIMERS__*/` placeholder. So the obvious way to show a new
primer design in the report is to re-run the pipeline -- which re-fetches
inDelphi and primer data for every one of the ~240 guides and, if the CHOPCHOP
job has expired, submits a new one and risks renumbering the guides (the whole
reason `refresh_job.py` exists). None of that has anything to do with primers.

The guide data in a generated report is already correct and does not depend on
the primer design at all. So this rewrites exactly one thing: the value of the
`const PRIMERS = ...;` literal, replaced with the current shortlist filtered to
that report's gene. Everything else in the file is left byte-identical.

The value's end is located by parsing it with json.JSONDecoder().raw_decode,
not by searching for a `;` or a newline -- the JSON blob is one long line and
contains both inside strings.

Usage
-----
  python inject_primers.py --check out_geneA/geneA_danRer11.html
  python inject_primers.py --apply out_geneA/geneA_danRer11.html out/yap1_danRer11.html
  python inject_primers.py --apply --candidates primers_out/candidates.json *.html
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys

MARK = "const PRIMERS = "


def find_blob(html: str) -> tuple[int, int, object]:
    """Locate the PRIMERS literal; return (start, end, parsed value)."""
    i = html.find(MARK)
    if i < 0:
        sys.exit("no `const PRIMERS = ` in this file -- is it a generated report?")
    if html.find(MARK, i + 1) >= 0:
        sys.exit("`const PRIMERS = ` occurs more than once -- refusing to guess")
    start = i + len(MARK)
    if html.startswith("null", start):
        return start, start + 4, None
    if html.startswith("/*__PRIMERS__*/", start):
        # the un-substituted template: nothing to replace, the pipeline does it
        sys.exit("this is the template, not a generated report "
                 "(the placeholder is still there)")
    try:
        val, end = json.JSONDecoder().raw_decode(html, start)
    except ValueError as e:
        sys.exit(f"could not parse the PRIMERS literal: {e}")
    if html[end:end + 1] != ";":
        sys.exit(f"the PRIMERS literal is not followed by ';' but by "
                 f"{html[end:end + 20]!r} -- refusing to edit")
    return start, end, val


def gene_of(html: str) -> str | None:
    """The gene this report was generated for, from the DATA payload."""
    m = re.search(r'"gene":\s*"([^"]+)"', html)
    return m.group(1) if m else None


def blob_for(cand: dict, gene: str) -> str:
    loci = [L for L in cand.get("loci", []) if L.get("gene") == gene]
    if not loci:
        sys.exit(f"candidates.json has no loci for gene {gene!r} "
                 f"(it has: {sorted({L.get('gene') for L in cand.get('loci', [])})})")
    return json.dumps({"shared_tm": cand.get("shared_tm"),
                       "params": cand.get("params"),
                       "tm_scan": cand.get("tm_scan", []),
                       "loci": loci}).replace("</", "<\\/")


def describe(val) -> str:
    if not val:
        return "none (null)"
    loci = val.get("loci", [])
    pairs = sum(len(L.get("pairs", [])) for L in loci)
    rec = sum(len(L.get("recommended", [])) for L in loci)
    p = val.get("params") or {}
    return (f"{len(loci)} loci, {pairs} candidates, {rec} recommended pairs, "
            f"Tm {val.get('shared_tm')} C, generated {p.get('generated', '?')}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reports", nargs="+", help="generated report HTML file(s)")
    ap.add_argument("--candidates", default=None,
                    help="candidates.json (default primers_out/candidates.json "
                         "next to this script)")
    ap.add_argument("--check", action="store_true",
                    help="say what would change, write nothing (default)")
    ap.add_argument("--apply", action="store_true", help="write the files")
    args = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    cpath = args.candidates or os.path.join(here, "primers_out", "candidates.json")
    with open(cpath, encoding="utf-8") as fh:
        cand = json.load(fh)
    print(f"candidates: {cpath}")
    print(f"            {describe(cand)}\n")

    changed = 0
    for path in args.reports:
        html = open(path, encoding="utf-8").read()
        gene = gene_of(html)
        if not gene:
            sys.exit(f"{path}: could not read the gene from the DATA payload")
        start, end, old = find_blob(html)
        new = blob_for(cand, gene)
        print(f"{path}  (gene {gene})")
        print(f"  now: {describe(old)}")
        print(f"  new: {describe(json.loads(new.replace(chr(92) + '/', '/')))}")
        if html[start:end] == new:
            print("  identical -- nothing to do")
            continue
        changed += 1
        out = html[:start] + new + html[end:]
        # Everything outside the literal must be untouched.
        assert out[:start] == html[:start] and out[len(out) - (len(html) - end):] == html[end:]
        if args.apply:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(out)
            print(f"  rewrote {len(html)} -> {len(out)} bytes "
                  f"(only the PRIMERS literal changed)")
        else:
            print(f"  would rewrite {len(html)} -> {len(out)} bytes")

    if not args.apply and changed:
        print(f"\n{changed} file(s) would change. Re-run with --apply.")
    elif args.apply:
        print(f"\n{changed} file(s) updated.")


if __name__ == "__main__":
    main()
