#!/usr/bin/env python3
"""
refresh_job.py -- re-point an existing report at a fresh CHOPCHOP job.

The problem
-----------
CHOPCHOP expires finished jobs after a while. Once it has, every
`.../results/<job>/details/<rank>` link in a generated report 404s, and the
"CHOPCHOP job" link in the header dies with them. The report's *data* is fine --
it was all captured into the HTML at generation time -- only the outbound links
are stale.

What this does, and does not do
-------------------------------
It submits a new job with the **identical** option array, waits for it, and then
matches the guides in the existing report to the guides in the new job by
genomic identity: (chrom, start, strand, protospacer). Not by rank. Rank is the
thing under suspicion, so it cannot also be the thing we trust.

Only if every guide matches AND every matched guide kept the same rank does it
rewrite anything, and then it rewrites exactly two strings in the embedded JSON:
`job` and `jobUrl`. Scores, SNP annotations, inDelphi numbers, primers, protein
context -- all untouched, because none of them are re-derived.

Starred guides live in localStorage under `gr-pins-<gene>-<build>`, which does
not contain the job id, so they survive this. But they are stored as *rank
numbers*, which is exactly why a rank mismatch has to abort rather than be
papered over: silently renumbering would move someone's starred set onto
different guides.

Usage
-----
  python refresh_job.py --check out_geneA/geneA_danRer11.html     # dry run, no writes
  python refresh_job.py --apply out_geneA/geneA_danRer11.html
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chopchop_pipeline import (  # noqa: E402
    CHOPCHOP, PipelineError, build_opts, get_text, parse_results, submit, wait_for,
)

PREFIX = "const DATA = "


def read_report(path: str) -> tuple[list[str], int, dict]:
    """Return (lines, index of the DATA line, parsed payload)."""
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if line.startswith(PREFIX):
            s = line.rstrip("\n")
            if not s.endswith(";"):
                sys.exit(f"{path}: DATA line does not end in ';'")
            return lines, i, json.loads(s[len(PREFIX):-1])
    sys.exit(f"{path}: no '{PREFIX}' line found")


def job_alive(job: str) -> bool:
    try:
        get_text(f"{CHOPCHOP}/results/{job}/results.tsv", timeout=30)
        return True
    except PipelineError:
        return False


def identity_old(g: dict) -> tuple:
    return (str(g.get("chrom")), int(g.get("start")), str(g.get("strand")),
            str(g.get("proto")).upper())


def identity_new(r: dict) -> tuple:
    """Build the same key from a raw results.tsv row.

    `Genomic location` is `chrom:start`. The pipeline stores that coordinate
    through verbatim into `start`, so no convention conversion belongs here --
    an earlier version of this function "helpfully" subtracted 1 and matched
    0/240. `Target sequence` is the 23-mer (protospacer + PAM); the report
    stores the 20-mer protospacer separately, so trim the PAM.
    """
    loc = r["Genomic location"]
    chrom, pos = loc.rsplit(":", 1)
    target = r["Target sequence"].upper()
    return (chrom, int(pos), r["Strand"], target[:20])


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Re-point a generated report at a fresh CHOPCHOP job.")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true",
                      help="submit, compare, report -- write nothing")
    mode.add_argument("--apply", action="store_true",
                      help="rewrite job/jobUrl if and only if the ranks agree")
    ap.add_argument("report", help="an existing *_danRer11.html")
    ap.add_argument("--job", default=None, metavar="ID",
                    help="compare against this already-finished job instead of "
                         "submitting a new one (re-checking should not cost "
                         "CHOPCHOP another full run)")
    ap.add_argument("--scoring", default="MORENO_MATEOS_2015")
    ap.add_argument("--pam", default="NGG")
    ap.add_argument("--product-min", type=int, default=150)
    ap.add_argument("--product-max", type=int, default=290)
    args = ap.parse_args()

    lines, idx, data = read_report(args.report)
    gene, build, old_job = data["gene"], data["build"], data["job"]
    guides = data["guides"]
    print(f"report   : {args.report}")
    print(f"gene     : {gene} / {build}")
    print(f"old job  : {old_job}")
    print(f"guides   : {len(guides)}")

    if job_alive(old_job):
        print("\nThe existing job is still alive on CHOPCHOP -- nothing to refresh.")
        return
    print("           (confirmed dead)")

    if args.job:
        new_job = args.job
        print(f"\nreusing the finished job {new_job} (no new submission)")
        tsv = get_text(f"{CHOPCHOP}/results/{new_job}/results.tsv", timeout=60)
    else:
        opts = build_opts(build, args.scoring, args.pam,
                          args.product_min, args.product_max)
        print(f"\nsubmitting a fresh job with the identical option array ...")
        new_job = submit(gene, opts)
        print(f"new job  : {new_job}")
        print("waiting for it to finish (this takes a few minutes) ...")
        tsv = wait_for(new_job)
    rows = parse_results(tsv)
    print(f"           new job returned {len(rows)} guides")

    new_by_id = {identity_new(r): int(r["Rank"]) for r in rows}
    matched, missing, moved = 0, [], []
    for g in guides:
        key = identity_old(g)
        if key not in new_by_id:
            missing.append(g["rank"])
            continue
        matched += 1
        if new_by_id[key] != int(g["rank"]):
            moved.append((g["rank"], new_by_id[key]))

    print(f"\nmatched by (chrom, start, strand, protospacer): "
          f"{matched}/{len(guides)}")
    print(f"missing from the new job                      : {len(missing)}")
    print(f"matched but renumbered                        : {len(moved)}")
    if missing:
        print(f"  missing ranks: {missing[:20]}"
              + (" ..." if len(missing) > 20 else ""))
    if moved:
        print(f"  moved (old -> new): {moved[:20]}"
              + (" ..." if len(moved) > 20 else ""))

    clean = (matched == len(guides) and not moved)
    if not clean:
        print("\nRANK MAPPING CHANGED -- refusing to rewrite.")
        print("Starred guides are stored as rank numbers, so renumbering here would")
        print("silently move your saved set onto different guides. The report keeps")
        print(f"working offline as-is; the new job {new_job} is live if you want to")
        print("browse it by hand.")
        sys.exit(2)

    print("\nEvery guide matched and kept its rank. Links can be re-pointed safely.")
    if args.check:
        print("(--check: nothing written)")
        return

    data["job"] = new_job
    data["jobUrl"] = f"{CHOPCHOP}/results/{new_job}/"
    # The report warns when its CHOPCHOP links are likely dead. That age has to
    # be measured from the refresh, not from `generated`, or a freshly
    # re-pointed report keeps claiming to be stale.
    data["job_refreshed"] = time.strftime("%Y-%m-%d %H:%M")
    blob = json.dumps(data).replace("</", "<\\/")
    lines[idx] = f"{PREFIX}{blob};\n"
    with open(args.report, "w", encoding="utf-8") as fh:
        fh.writelines(lines)
    print(f"rewrote {args.report}: job -> {new_job}")


if __name__ == "__main__":
    main()
