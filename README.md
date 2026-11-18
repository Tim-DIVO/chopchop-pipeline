# chopchop_pipeline

End-to-end CRISPR guide selection for zebrafish F0 knockouts, following
Kroll & Rihel's protocols. Submits (or reuses) a CHOPCHOP job, pulls everything
CHOPCHOP computes, annotates it against the RefSeq gene model, checks strain SNPs,
and writes an annotated CSV plus an interactive HTML report.

## What is in this folder

```
chopchop_pipeline.py    guide selection: the main script
protein_annot.py        REQUIRED - UniProt domains, CDS geometry, in-frame ATGs
protein_align.py        REQUIRED - BLOSUM62 affine-gap global alignment
report_template.html    REQUIRED - must sit next to the script
transcript_check.py     is a frameshift in exon 1 escapable via another isoform?

primer_design.py        genotyping primers: one shared Tm across the plate
verify_primers.py       independent structural check of what it produced
ispcr_scan.py           genome-wide specificity from a local .2bit
blast_specificity.py    the same question via NCBI, when a genome is not to hand
inject_primers.py       put a new primer design into an existing report
refresh_job.py          re-point a report at a live CHOPCHOP job

README.md               this file
```

The guide half and the primer half are independent: `chopchop_pipeline.py`
decides *where to cut*, everything in the second block decides *how to read the
cut back out*. They meet only at `candidates.json`, which the report renders and
`inject_primers.py` can swap in without re-running the pipeline.

## Install

Nothing to install. Python 3.9+, standard library only.

## Use

```bash
# from scratch -- submits the job for you
python chopchop_pipeline.py --gene yap1 --build danRer11 -o out

# reuse a job you already ran in the browser
python chopchop_pipeline.py --job-url https://chopchop.cbu.uib.no/results/<id>/ -o out
```

Outputs `out/<gene>_<build>.csv` and `out/<gene>_<build>.html`.

### Options

| flag | default | meaning |
|---|---|---|
| `--gene` | — | gene symbol or Ensembl ID. Needed unless `--job-url` is given |
| `--job-url` | — | reuse an existing CHOPCHOP results URL instead of submitting |
| `--build` | `danRer11` | UCSC assembly |
| `--scoring` | `MORENO_MATEOS_2015` | efficiency model — the one Kroll recommends for zebrafish |
| `--snp-top` | `0` (all) | query SNPfisher for the top N by CHOPCHOP rank; `0` = all |
| `--product-min/max` | `150` / `290` | primer3 amplicon size window |
| `--preselect` | — | guides to pre-save and pin to the top of the report: CHOPCHOP ranks and/or pasted `.../details/<n>` URLs, comma- or whitespace-separated |
| `--delay` | `0.2` | seconds between per-guide fetches |
| `-o, --out` | `chopchop_out` | output directory |

`--snp-top` defaults to **all guides**, because SNPfisher is queried in a **single
batched request** regardless of how many windows — restricting it saves essentially
nothing. Set a number only if you want to keep the request small. Every guide gets its
Zv9 coordinates computed either way, so the HTML's live SNP button works on any row.

## What the columns mean

| column | source |
|---|---|
| `efficiency` | CRISPRscan (Moreno-Mateos 2015) — *will Cas9 cut here* |
| `frameshift` | inDelphi (Shen 2018) — *will the repair shift the frame*. 66% = random, below that is worse than chance; Kroll wants >80 |
| `precision`, `mh_del`, `exp_indel` | further inDelphi repair statistics |
| `exon`, `symmetry`, `bp_to_splice` | computed here from the UCSC RefSeq model, using the longest-CDS transcript as reference, per Kroll |
| `mm0`–`mm3` | CHOPCHOP off-targets, **ungapped mismatches only** |
| `zv9`, `snps` | liftOver to Zv9 + SNPfisher strain SNPs across the 23-mer |
| `snp_fli`, `snp_tl`, `snp_wik` | of those, how many are actually segregating in each line. SNPfisher genotyped only these three; for AB or TU use the `snps` total |
| `snp_pam`, `snp_seed`, `snp_distal` | of the SNPs in the 23-mer, how many sit on a **G of the NGG PAM**, in the **PAM-proximal seed** (protospacer 11-20), or in the tolerated remainder (protospacer 1-10 plus the N of NGG) |
| `snp_worst`, `snp_worst_zone` | the severity of the worst-placed SNP: `critical` (PAM G), `high` (seed), `low` (distal), `negligible` (the N) |
| `snp_pam_fli/tl/wik` | PAM-G SNPs actually segregating in each genotyped line |
| `aa`, `aa_total` | the residue the cut site falls in, and the protein length |
| `domains` | annotated domains spanning that residue. Entries marked `(projected)` were mapped from the human ortholog through the alignment rather than annotated on the zebrafish protein |
| `cons`, `cons_sim` | percent identity and BLOSUM-positive similarity to the human ortholog across the 11 residues centred on the cut |
| `next_atg`, `atg_dist`, `reinit_keeps` | the next **in-frame ATG** downstream, how many codons away it is, and how much of the protein a re-initiated product would keep |
| `first_coding_exon` | which exon carries the start codon, so the report can mark cuts in it |
| `primer1..5_*` | CHOPCHOP/primer3 pairs: sequences, Tm, product size |

### Two traps this encodes

**MM0 = 1 on every guide is usually an assembly ghost copy.** Some loci appear twice
in the reference — once on the real chromosome, once on an alternate contig
(`ALT_CTG18_1_40`, `chr9_KZ115297v1_alt`, …). That is a bookkeeping artefact, not a
second cut site. Judge on MM1–MM3. The HTML's off-target filter defaults to
"MM1–MM3 all zero" for this reason.

**The off-target CSV you can download from CHOPCHOP is not the same search as the
MM columns.** It permits gapped alignments (bulges), so it lists far more hits. A
"distance 2" row there may be one mismatch plus a 1-nt deletion, which is why a guide
can show MM2 = 0 and still appear in that file.

## HTML report

Self-contained, works offline, opens in any browser.

- sortable columns, live filters (efficiency, frameshift, splice distance, SNPs,
  self-complementarity, off-targets, symmetry, exon chips, sequence search)
- **adjustable ranking weights** — each factor is min–max normalised across the
  guides currently shown, weighted, and summed. Set a weight to 0 to ignore it.
- **a SNP-position penalty that is multiplicative, not additive.** A SNP on a G of
  the NGG PAM does not make a guide slightly worse, it means the guide probably
  does not cut, so the score is multiplied by 0.05 at the default weight rather
  than docked a few points. A seed SNP costs about 29%. Both follow the strain
  selector, so a variant absent from your line is not counted.
- **protein context per guide** — the residue hit, the domains spanning it, and
  percent identity to the human ortholog at that point, with a *Protein context*
  card listing every domain span and the overall identity
- **alternative-start-codon awareness** — the next in-frame ATG downstream of each
  cut, how far away it is, and how much protein a re-initiated product would keep,
  plus a filter for "skip the first coding exon"
- **saved guides** — click the ☆ on any row to pin it to the top of the table. The
  list is kept per gene in the browser, survives a reload and a re-run, and can be
  exported in bulk: CHOPCHOP links, bare 20-mers, IDT-ready FASTA, a CSV, or a
  links `.txt`. `--preselect` seeds the list from the command line.
- click any row for primers, full inDelphi stats, position, protein and SNP detail
- **CC** opens that guide on CHOPCHOP; **SNP** re-queries SNPfisher live in a new tab
  (a real form POST, so it works for rows the script did not pre-query); **⧉** copies
  the 20 nt without the PAM, ready for IDT
- a **wild-type strain** selector (any / TL / WIK / fli1a-eGFP) that switches the
  SNP column, filter, scoring penalty and summary to that line's allele
  frequencies. `0 /2` means zero SNPs that matter in your strain, two in total
- export the filtered set back out to CSV
- a second **Primers** tab: the shortlist per locus with the ordered pairs
  pre-ticked, each pair's predicted in-silico PCR result, a conditions panel
  recording what the search ran under, and order sheets with or without the
  Illumina overhangs

## How it works — the full logic

Every number in the output is either fetched from a named source or derived by the
rules below. Nothing is estimated or inferred.

### 1. The CHOPCHOP job

Submitted as `POST /` with a JSON body `{opts, geneInput, fastaInput, isIsoform,
forSelect}`, which returns `{jobId}`. The `opts` array mirrors what the web form
sends. Three entries matter:

- `-scoringMethod MORENO_MATEOS_2015` — CRISPRscan, the efficiency model Kroll
  recommends for zebrafish, and the only one benchmarked to correlate with in-vivo
  results (Uribe-Salazar 2022, and even then only r = 0.27–0.31).
- `-repairPredictions mESC` — without this, no inDelphi data is generated at all.
- `-P` with `-3 PRODUCT_SIZE_MIN=…` — without these, no primers.

Identical queries return the **same jobId**, so re-running a gene costs nothing.
The script then polls `results.tsv`, which 404s until the job finishes.

### 2. Verifying every guide against the genome

CHOPCHOP reports a 23-mer and a 1-based start coordinate, but does not say whether
that 23-mer reads forward or reverse at that position. The script fetches the locus
from the UCSC API and matches each guide both ways:

```
window == guide          -> plus strand
revcomp(window) == guide -> minus strand
neither                  -> hard error (wrong assembly)
```

This is a genuine check, not a formality: it is what catches a job run on a
different genome build.

### 3. Cut site

Cas9 makes a blunt cut 3 bp from the PAM-proximal end of the protospacer. The script
records the base immediately left of the break:

```
plus  strand:  cut = start + 16
minus strand:  cut = start + 5
```

Validated independently: CHOPCHOP's inDelphi output contains a 120-nt reference
sequence centred on its own cut site, and the script confirms it matches the genome
at the computed position. (An earlier off-by-one on the minus strand was found
exactly this way.)

### 4. Reference transcript

Kroll: *"Choose the longest protein-coding transcript as reference for numbering the
exons."* Longest **protein**, so transcripts are ranked by true coding length —

```python
sum(max(0, min(exon_end, cds_end) - max(exon_start, cds_start)) for each exon)
```

— **not** by `cdsEnd - cdsStart`, which is the genomic span and includes introns.
Those give different answers; ranking by genomic span is effectively arbitrary when
two transcripts share CDS boundaries.

Ties are common, because a curated `NM_` and a predicted `XM_` often describe the
same protein. The tiebreak prefers the **curated** record, then falls back to the
transcript name so the choice is deterministic.

### 5. Exon symmetry

An exon and the coding sequence are two independent intervals on the chromosome. The
coding part of an exon is their overlap:

```python
coding = max(0, min(exon_end, cds_end) - max(exon_start, cds_start))
```

Then:

| coding length | verdict |
|---|---|
| `% 3 != 0` | **asymmetric** — preferred; skipping the exon also shifts the frame |
| `% 3 == 0` | **symmetric** — skipping it deletes whole codons, frame survives |
| `== 0` | **noncoding** — the exon is pure UTR |

Using the **coding** length rather than the whole exon is what makes this correct, and
Kroll says so explicitly: *"the first and last exon may include a UTR… do not count
the length of the UTR."* Only the first and last exons are ever clipped, because
middle exons sit entirely inside the CDS.

It changes answers. For zebrafish *yap1*, exon 1 is 603 bp (not divisible by 3, so
naïvely asymmetric) but only 201 bp of that is coding (divisible by 3, so actually
symmetric). Exon 8 flips the same way: 2448 bp total, 239 bp coding. Two of eight
exons get the wrong verdict if the UTR is not removed.

Cross-checked against Ensembl: `GET /overlap/id/<transcript>?feature=cds` returns the
per-exon CDS segments directly, and they match these computed values exactly.

### 6. Exon vs intron

The exon assignment uses the **cut site**, not the protospacer. A guide can overlap an
exon while cutting a few bp past the splice site — an indel there lands in the intron
and usually will not shift the reading frame, which is what Kroll warns about when he
says not to pick guides *"very close to the intron boundary"*.

Those guides get `region = intron` and a **negative** `bp_to_splice` (how far past the
boundary the break falls), so any filter of `bp_to_splice >= 0` excludes them. For
*yap1*, 22 of 241 guides are in this category.

Kroll gives **no numeric margin** for "very close". The default filter of 0 only
excludes genuinely intronic cuts; any positive threshold is your own judgement call.

### 7. Off-targets

`MM0`–`MM3` come straight from CHOPCHOP and count **ungapped mismatches only**.

Two things are easy to misread:

- **`MM0 = 1` on every guide is usually an assembly ghost copy.** Some loci appear
  twice in the reference — once on the real chromosome, once on an alternate contig
  (`ALT_CTG18_1_40`, `chr9_KZ115297v1_alt`, …). It is a bookkeeping artefact, not a
  second cut site. Kroll's test is a suspicious contig name that recurs across many
  guides. Judge on MM1–MM3; the HTML's off-target filter defaults to that.
- **CHOPCHOP's downloadable off-target CSV is a different, more permissive search**
  that allows gapped alignments (bulges). A "distance 2" row there can be one mismatch
  plus a 1-nt deletion, which is why a guide can read `MM2 = 0` and still appear in it.

### 8. inDelphi

Fetched per guide from `<rank>_repStats.json`. inDelphi (Shen 2018) predicts the
*distribution of repair outcomes* at a cut, driven largely by microhomology in the
flanking sequence. `Frameshift frequency` is the share of outcomes whose length is not
a multiple of 3.

Calibration: **66 % is what random indel lengths would give**, so anything below that
is worse than chance. Kroll wants >80 %. Note the model was trained on mouse embryonic
stem cells, not zebrafish embryos.

Efficiency and frameshift answer different questions — *will it cut* versus *will the
repair shift the frame* — and a guide can score well on one and badly on the other.

### 9. Liftover to Zv9

SNPfisher was built in 2015 on Zv9 (= `danRer7`); CHOPCHOP reports GRCz11
(= `danRer11`). The script converts with the UCSC chain file
`danRer11ToDanRer7.over.chain.gz`, parsed directly — no `liftOver` binary, no
`pyliftover` dependency.

Every mapping is **validated by sequence**: the 23-mer is fetched from both builds and
compared. A mismatch is reported as `liftover-seq-mismatch` and skipped rather than
queried. This is strictly better than the BLAST-and-eyeball step in the protocol,
because it is automatic and checkable.

### 10. SNPfisher

One batched POST to `snpTextOracle.jsp` with all windows at once (the form accepts one
locus per line), so the number of guides queried barely affects cost.

**The trap:** the form ships its three frequency comparators as `= 0.0`, which matches
only SNPs at exactly zero frequency in *all three* strains — i.e. it silently returns
nothing. The script forces them to `>=`.

Results are split per line (`snp_fli`, `snp_tl`, `snp_wik`) as well as a total,
because whether a SNP matters depends on which wild-type strain you actually use.
SNPfisher genotyped only Tg(fli1a-eGFP), WIK and TL — for AB or TU, use the total and
treat it as weaker evidence.

### 11. Ranking in the HTML

The `Score` column is **not** from CHOPCHOP. Each factor is min–max normalised across
the guides currently visible, multiplied by its weight, summed, and divided by the
total weight. Penalties (off-targets, SNPs, self-complementarity) enter as `1 - x`.
Setting a weight to 0 removes that factor. Defaults weight frameshift highest (1.5),
then efficiency and off-targets (1.0).

Two factors are applied **multiplicatively** instead, after that sum: a SNP on a G
of the PAM and a SNP in the PAM-proximal seed. See section 12.

This is a convenience for exploring trade-offs, not a published metric. The underlying
CHOPCHOP rank is kept in its own column.

### 12. Where a SNP sits in the 23-mer

The original SNP column counts variants anywhere across the target. That treats a
variant on the PAM and one at the far 5' end of the protospacer as the same problem,
and they are not.

**How the position is worked out.** The Zv9 window stored per guide always runs low
to high in genomic coordinates, while `target` is written 5'→3' on the *guide*
strand. So for a plus-strand guide the k-th base of the target sits at
`window_start + k`, and for a minus-strand guide at `window_end - k` — which puts
the PAM at the three **lowest** coordinates. Each SNPfisher hit is converted to that
offset `k` and binned:

| offset in the 23-mer | zone | severity | why |
|---|---|---|---|
| 22, 21 | the two **G**s of NGG | **critical** | Cas9 reads these two bases directly, through Arg1333 and Arg1335 (Anders 2014). A variant here normally costs the cut. If the change happens to make NAG, roughly a fifth of the activity can survive (Hsu 2013; Jiang 2013) — still not a guide you want. |
| 20 | the **N** of NGG | negligible | not read by the protein |
| 10–19 | protospacer 11–20, the **seed** | high | single mismatches in the ~10 nt nearest the PAM usually leave well under a quarter of the activity (Hsu 2013; Doench 2016) |
| 0–9 | protospacer 1–10, PAM-distal | low | single mismatches here are largely tolerated |

The penalty is multiplicative because the failure is not graded: `score × 0.05^(w/3)`
for a PAM-G hit and `× 0.5^(w/3)` for a seed hit, so at the default weights (3 and
1.5) a PAM-G hit costs 95% of the score and a seed hit 29%. Setting either weight to
0 switches that penalty off. Both follow the strain selector: for a named line only
variants segregating in that line count, because a variant fixed in WIK says nothing
about a TL fish. On *any strain* the raw counts are used, which is the conservative
reading.

SNPfisher reports per-strain frequencies but not a per-strain breakdown by position,
so the seed count for a named strain is gated on that strain carrying at least one
variant in the window rather than resolved variant by variant. The per-variant
frequencies are all in the expanded row.

### 13. Protein context: domains and conservation

Kroll's protocol deliberately avoids domain-aware design, because it is meant for
screens of poorly annotated genes. For two named genes the annotation exists, so the
report shows it rather than ignoring it.

1. **The CDS is spliced from UCSC**, exon by exon, in transcription order. On the
   minus strand each block is reverse-complemented *individually* and the blocks are
   then joined. Reverse-complementing the finished concatenation instead would also
   reverse the order of the exons and yield a sequence of exactly the right length
   that does not begin with ATG — a silent failure, so the pipeline checks that the
   spliced CDS starts with ATG and is a multiple of 3, and skips the analysis loudly
   if not.
2. **Each cut site is converted to a residue number** through the per-exon coding
   offsets, so a guide reports the codon it actually breaks, not just its exon.
3. **Domains come from UniProt** for the zebrafish protein, plus the human
   ortholog's features projected onto zebrafish residue numbers. Projected spans are
   labelled as such everywhere; their edges are only as good as the alignment.
4. **The alignment is computed locally** — Needleman–Wunsch with affine gaps (Gotoh)
   on BLOSUM62, open −11, extend −1, the BLAST and Clustal protein defaults. For two
   sequences this is the same computation Clustal Omega reduces to, so there is no
   external service to call and no e-mail address to hand over. It lives in
   `protein_align.py` and is standard library only, like the rest.
5. **Conservation per guide** is percent identity, and BLOSUM-positive similarity,
   over the 11 alignment columns centred on the cut residue.

What the **annotated-domain bonus** in the weights panel does and does not mean: it
adds a small, adjustable bonus when the cut falls inside a domain whose name matches
a short list of functional keywords (Tudor, tether, TEAD, transactivation, reductase,
WW, homeobox, kinase, DNA binding, zinc finger). For a frameshift the distinction
between cutting *inside* a domain and cutting *before* it is thin — everything
downstream is lost either way — so treat this as a weak "is this part of the protein
annotated as doing something" signal, not as domain-targeted design in the
Hoshijima 2019 sense. The residue number and the domain column are the parts to read;
the bonus is only there so the sort can be nudged. Its default weight is 0.3.

A few domains that matter for design are not in UniProt — the YAP1 TEAD-binding
domain, its PDZ-binding motif, the LBR sterol-reductase domain. Those are curated in
`EXTRA` in `protein_annot.py` in **human** coordinates, each with its source and a
one-line note on why it changes a design decision, and are projected like any other
human feature.

### 14. Alternative start codons

The last item on Kroll's exon-priority list is "not exon 1", because a 5' frameshift
can be bypassed: the ribosome re-initiates at a downstream in-frame ATG and makes an
N-terminally shortened protein. This is documented in stable zebrafish lines
(Anderson 2017; Lalonde 2017) and in about a third of 193 HAP1 knockout lines with
residual protein (Smits 2019).

The pipeline makes that concrete rather than leaving it as a rule of thumb. It scans
the spliced CDS for every in-frame ATG after the initiator and, for each guide,
reports the next one downstream, its distance in codons, and the fraction of the
protein a product starting there would keep. The re-initiation penalty decays with
distance and is gone past 150 codons, because a cut far upstream of the nearest
alternative start has no cheap escape route. The exon chips name the in-frame ATGs
inside each exon, so the first coding exon can be judged rather than simply avoided.

### Why UCSC rather than Ensembl

Both serve the same assembly (`danRer11` = `GRCz11`) and, for *yap1*, the same gene
model — RefSeq `NM_001139480.1` and Ensembl `ENSDART00000098914` both give 8 exons and
442 aa. Ensembl can even return per-exon CDS segments directly, skipping the
intersection in step 5.

UCSC is used because everything else in the pipeline is already in UCSC coordinates:
CHOPCHOP's output (`chr18:37351408`), the liftOver chain, and the sequence API.
Ensembl names the same chromosome `18` rather than `chr18`, so mixing the two would
add a conversion layer for no gain.

The caveat: RefSeq and Ensembl are **independently curated** and do differ for some
genes. Kroll's protocol numbers exons from Ensembl, this pipeline numbers them from
RefSeq, and for a gene with genuinely divergent isoforms the two can disagree. The
report names the transcript it used; if exon numbers matter for a specific gene, check
it against Ensembl's transcript table.

## Caveat: undocumented endpoints

The per-guide endpoints (`<rank>_repStats.json`, `primer_<rank>.json`) and the job
submission payload are **not a supported public API**. They were read out of
CHOPCHOP's own client-side JavaScript (`detailsController.js`, `formController.js`).
They work as of 2026-11 but can break without notice.

The script is written to fail loudly rather than emit blank columns — for instance it
raises if primer files are fetched but no pairs parse, which is exactly how a silent
layout change would otherwise slip through.

Similarly, SNPfisher's frequency comparators ship as `= 0.0`, which matches only SNPs
at exactly zero frequency in all three strains and so silently returns nothing. The
script forces them to `>=`.

## Sources

- Kroll, *How to select the best gRNA(s) for frameshift knockouts in zebrafish*,
  dx.doi.org/10.17504/protocols.io.81wgb6r5qlpk/v1
- Kroll & Rihel, *F0 knockout — single gene* v3,
  dx.doi.org/10.17504/protocols.io.5qpvo52wdl4o/v3
- Kroll et al. 2021, eLife 10:e59683
- Moreno-Mateos et al. 2015, Nat Methods (CRISPRscan)
- Shen et al. 2018, Nature (inDelphi)
- LaFave et al. 2014, Genetics (SNPfisher)

## Cache

The liftOver chain (~2 MB) and the SNPfisher form are cached in `.cache/` next to
the script, so they download once per machine rather than once per gene. Safe to
delete; it will be refetched.

---

# Companion tools

The guide pipeline above picks *where to cut*. These two answer the next two
questions. Both are deterministic; the judgement that goes with them lives in
the `crispr-primers` skill (`.claude/skills/crispr-primers/SKILL.md`), which is
deliberately kept separate so the scripts stay reproducible.

## primer_design.py -- genotyping primers for MiSeq

CHOPCHOP's Primer3 pairs do not make an amplicon panel. In order of how much
each actually mattered, measured on a twelve-guide test set across two genes:

1. **Length.** Only ~32% of the pairs CHOPCHOP returned were <= 233 bp
   (401/1200 and 374/1205 on the two genes), and for six of the twelve guides all
   five returned pairs were too long. This is the one that blocked the work.
   Note that re-running CHOPCHOP with `--product-max 233` would have addressed
   this alone -- the rest is what a re-run would not have fixed.
2. **Known variants are not excluded.** A SNP under a primer's 3' end stops
   extension, and the PCR then fails on some fish and looks like a failed
   injection. ~21% of zebrafish target sites carry an annotated SNP.
3. **Cut position inside the product is not constrained.** CHOPCHOP only
   requires the guide to be somewhere in the amplicon; CRISPResso2 discards the
   terminal 15 bp of every read by default.
4. **Shared Tm across loci** -- a convenience, and a smaller gain than it might
   look. Measured spread of CHOPCHOP's Tm values was 57.2-62.1 C on one gene
   and 58.5-61.5 C on the other, which one gradient PCR would have settled.

An earlier version of this file said left and right primer "can differ by 6 C",
citing CRISPOR's `runPrimer3()` leaving `PRIMER_PAIR_MAX_DIFF_TM` at Primer3's
default of 100. That is true of the *setting* but misleading about the
*behaviour*: Primer3 optimises toward `PRIMER_OPT_TM`, so pairs come out close
anyway. In CHOPCHOP's actual output here the median left/right difference was
**0.30 C**, with only 6% of pairs above 1 C. Do not repeat the 6 C claim.

This script treats the whole guide set as one design problem.

```bash
pip install primer3-py certifi

python primer_design.py \
  --from out_geneA/geneA_danRer11.csv:3,6,34,15,33,5 \
  --from out/yap1_danRer11.csv:5,16,17,13,10,21 \
  -o primers_out
```

Pass **every locus in one invocation** -- the shared Tm is optimised across the
whole set, so splitting the genes defeats the point.

What it does, in order: fetches the template from UCSC (cached on disk, with
backoff -- the REST API returns 429 if you re-run without it); collects variants
from Ensembl (GRCz11) and SNPfisher (via the Zv9 liftOver chain, lifting every
base individually rather than interpolating); scans Tm in 0.25 C steps for the
value that satisfies the most loci; then designs with Primer3 and walks a fixed
relaxation ladder, reporting which rung each locus needed.

Constraints, and where they come from:

| constraint | default | provenance |
|---|---|---|
| product length | 140-233 bp | the MiSeq budget; Kroll et al. 2021's own panel was 151-198 bp across 42 loci |
| cut to product end | >= 40 bp | Kroll's panel ran 39-93 bp; CRISPResso2 ignores indels in the terminal 15 bp by default |
| left/right Tm difference | <= 0.5 C | the thing CRISPOR leaves unset |
| shared Tm across all loci | scanned | so one annealing temperature runs the plate |
| variants under a primer | none | Ensembl + SNPfisher; a variant in the 3' 8 nt is never tolerated at any tier |
| primers on the protospacer | forbidden | `SEQUENCE_TARGET` spans guide + flank |
| primer length | 18-27 nt, Primer3 aims for 23 | length is the cheapest specificity there is, measured -- see below |

**Why the oligos are long.** An exact 20-mer is expected to recur in a 1.4 Gb
genome about 4^3 = 64 times more often than an exact 23-mer, and a longer oligo
costs nothing but synthesis. So the shortlist is ordered to prefer length
(`--sort length`, the default) and Primer3 is told to aim near the top of the
range. The limit of that argument, stated plainly: the Tm window is fixed at one
value for the whole plate, so a longer primer is necessarily more AT-rich, and
an AT-rich primer tolerates a *mismatched* partial site about as well as a short
GC-rich one. Length buys uniqueness of the exact match, not thermodynamic
discrimination. It is a good default, not a replacement for a specificity check.

The settings are the measured optimum on a twelve-guide test set, not a guess
(sweep run 2026-11-16, Tm scan + full design per configuration):

| min/opt/max | shared Tm | mean oligo | tiers | independent backups |
|---|---|---|---|---|
| 18/20/27 (old) | 57.75 C | 20.6 nt | 11 strict + 1 flank-30 | 6/12 |
| **18/23/27** | **59.75 C** | **22.9 nt** | **11 strict + 1 flank-30** | **7/12** |
| 18/25/27 | 59.75 C | 23.0 nt | 11 strict + 1 flank-30 | 5/12 |
| 20/25/27 | 59.75 C | 23.2 nt | 10 strict, hardest locus down a rung | 3/12 |
| 22/26/30 | 61.50 C | 24.6 nt | 9 strict, **hardest locus fails entirely** | 1/12 |

Past opt 23 the mean length stops improving while the candidate pool thins -- a
longer primer inside a fixed Tm window must be AT-rich, and there are fewer
places to put one. Raising `min` above 18 is the real trap: the hardest locus
in that set needs one short GC-rich oligo and has no pair at all without it. The higher shared Tm that
the longer primers bring (59.75 C) is itself a specificity gain: a more
stringent anneal.

**Shortlist order.** Rank 1 is chosen by: dTm bucketed to 0.1 C, then variants
under either primer, then the shorter oligo's length (to 2 nt), then product
length, then Primer3's penalty. Two quantities are bucketed on purpose -- a
0.05 C dTm difference predicts nothing about a PCR on a single block, and one
base of primer length should not outrank 59 bp of amplicon, so neither is
allowed to dominate a real difference. Variants rank above all of it: a SNP
under a primer risks allele dropout in exactly the genotyping assay these
primers exist for. (`--sort penalty` restores the old dTm-then-Primer3 order.)

Two failure modes it refuses to paper over:

- **A failed variant lookup is not "no variants".** If Ensembl errors out after
  retries, the locus is reported as `UNVERIFIED-<tier>` rather than clean.
- **Specificity is off by default and says so.** UCSC's web In-Silico PCR sits
  behind a Cloudflare bot check that explicitly points automated users at their
  FAQ; getting past it means spoofing a User-Agent, so the tool does not. Use
  `--ispcr-bin`/`--genome` with a local UCSC `isPcr`, run `blast_specificity.py`
  (below), or paste the pairs in `specificity_todo.txt` into NCBI Primer-BLAST.

### What `primers.csv` contains

`primers.csv` carries **one row per ordered pair**, primary and backup, not one
per locus: the file people order from must contain every oligo being ordered.
`id` (`<locus>_p<rank>`) is the unique key, since `name` now repeats, and the
downstream checkers key on it. `role`, `independent` and `backup_why` travel on
both rows of a locus, because a row gets copied out of its spreadsheet and a
non-independent backup must not travel without the sentence saying so.

The Illumina two-step overhangs (33 nt + 34 nt = 67, from doc 15044223 Rev. B)
are *not* added to the ordered sequence by default -- `primers.csv` carries them
in the `*_with_overhang` columns, and `--overhang-in-order-sheet` fuses them.
They do not consume sequencing cycles: the read starts after the overhang.

## transcript_check.py -- is exon 1 safe to target?

The main pipeline's `next_atg` answers the *within*-transcript escape route:
could translation re-initiate at a downstream in-frame ATG. This answers the
*between*-transcript one: if the gene has a coding isoform starting at a
different promoter, a frameshift in the reference transcript's first coding exon
is simply absent from that isoform's mRNA.

```bash
python transcript_check.py --gene yap1 --chrom chr18 \
       --from out/yap1_danRer11.csv:5,16,17,13,10,21
```

It prints every coding isoform with its CDS length, start and stop codon, and
exon count, says whether they share a start codon, and then reports for each
guide how many isoforms its cut site actually lands in.

On both test genes there are two annotated coding isoforms sharing one start
codon, so exon 1 is not disqualified by isoform structure in either. It did
catch one thing worth knowing: **yap1 rank 13 cuts in an exon that
`XM_005173662.4` skips**, so an indel there would be spliced out of that
transcript.

## verify_primers.py -- independent structural check of the primers

`primer_design.py` writes a CSV asserting where each primer sits, which guide it
flanks, and how long the product is. Those assertions come from the program that
made the choices, so a coordinate bug would produce a CSV that is wrong and
internally consistent. This script re-derives all of it and shares nothing with
the designer:

| | designer | verifier |
|---|---|---|
| sequence source | UCSC REST | **Ensembl REST** |
| coordinate convention | 0-based half-open | **1-based inclusive** |
| Tm | Primer3 2.6.1 | **SantaLucia 1998 written out longhand** |
| shared code | -- | **none** (no import from `primer_design`) |

```bash
python verify_primers.py primers_out/primers.csv \
    --guides geneA=out_geneA/geneA_danRer11.csv \
    --guides yap1=out/yap1_danRer11.csv \
    -o primers_out/verification.csv
```

Per pair it checks that the forward primer occurs **exactly once** in a window
around the cut, that the reverse primer's reverse complement occurs exactly once
downstream of it, that the measured product length equals the claimed one, that
the cut falls inside the product at the claimed offsets, and that neither primer
touches the protospacer. Exit code is 0 only if every locus passes.

Current result: **24/24 pass** (12 primary + 12 backup pairs), products
identical to the claimed values, cut offsets agreeing within the 1 bp allowed
for the half-open/closed convention.

Two caveats worth keeping in mind when reading a pass:

- The independent Tm agreeing to within 0.05 C is a **weak** check, not a strong
  one: the same nearest-neighbour table and the same salt correction were chosen
  deliberately, so it confirms the arithmetic, not the model. The sequence and
  position checks are the ones carrying weight.
- Specificity is **not** checked -- the window is local by construction. That
  still needs Primer-BLAST or a local `isPcr`.

The verifier was validated against a sabotaged input (one substituted base, one
inflated product length, one swapped primer pair, one cut moved 300 bp): all
four were caught with the correct diagnosis while the untouched rows still
passed. Re-run that negative control if the script is ever modified -- a checker
that only ever says PASS is worthless.

## ispcr_scan.py -- genome-wide specificity, locally

Primers that amplify two places genotype neither. Every route to checking that
was blocked, and the reasons are worth keeping because they will come back:

| route | why it is not used |
|---|---|
| UCSC web In-Silico PCR (`hgPcr`) | Cloudflare bot check that points automated clients at the downloads FAQ; getting past it means lying about the User-Agent |
| UCSC `isPcr` binary | distributed for Linux/macOS; this workstation is Windows |
| NCBI Primer-BLAST | no API, no batch mode |
| NCBI BLAST URL API | works, but the public queue is the bottleneck -- see below |

**What the BLAST slowness actually was.** The 2026-11-13 run was blamed on the
database being too large. It was not. Measured 2026-11-16: a **two-primer**
search submitted to the zebrafish-only genome database (`GPIPE/7955/current/
all_top_level`, word size 11, E 10) and an identical one against
`refseq_representative_genomes` were submitted within seconds of each other.
NCBI's own RTOE estimates were 5 s and 12 s. Both were still `WAITING`
**25 minutes later**, when the probe gave up. The queue, not the search, is the
cost, and it is not something a parameter fixes.

Running it on this design confirmed the measurement rather than the primers:
submitted 10:30, **still `WAITING` when it gave up at its 90-minute ceiling**,
having never left the queue. The RID is saved next to the output, so nothing is
lost -- but the local scan had answered the same question for all 332 candidate
pairs in the meantime, twice.

Three things were fixed anyway, because they were real and cheap:

- **A transport error no longer discards the wait.** The run that died at 32
  minutes died of `getaddrinfo failed` -- one DNS hiccup inside the poll loop.
  Every call now retries with backoff, and the RID is written to disk before the
  first poll and reused on the next run if the query set is unchanged, so a
  queued search survives a crash, a timeout, or a closed laptop.
- **Polling moved from 30 s to 60 s**, which is the floor NCBI asks for. A
  client that ignores the usage guidelines is one NCBI is entitled to throttle,
  and that is a candidate explanation for queue times far past the estimate.
- **The default database is now the zebrafish genome**, which is the right
  search space regardless of whether it is the bottleneck.
- **Query names are unique** (`<locus>_p<rank>`, not the locus name). BLAST
  groups its output by query-def, so the two pairs of a locus sharing a name
  merged their hits -- and pairing a primary forward primer with a backup
  reverse primer then invents products nobody ordered. This was introduced and
  caught in the same session by a patch applied without checking that it
  matched; the file ended up with the FASTA built from one key and the results
  read back under another, which would have reported *no hits at all* for every
  pair. **Assert the match count in a string-replacement patch.**

### What the local scan does

`ispcr_scan.py` reads a UCSC `.2bit` directly and reimplements isPcr's model of
priming, because that model is the right one:

1. **Perfect 3' seed** -- the last 15 bases (`--min-perfect`, isPcr's default)
   must match exactly. A mismatch in the last few bases stops extension; one at
   the 5' end usually does not. This is the filter that separates "looks a bit
   like the primer" from "will amplify", and counting raw BLAST hits without it
   over-reports off-targets wildly.
2. **Tolerant 5' tail** -- up to `--max-mismatch` (default 3) mismatches 5' of
   the seed. No gaps: a bulged primer that still primes is rare enough that
   modelling it would add more false positives than it removes.
3. **Pair the sites up** -- a plus-strand site and a minus-strand site on the
   same sequence, facing each other, within `--max-product` (default 3000 bp).
   Either primer may play either role at an off-target site, so both
   assignments are tried.

Soft-masked repeat blocks are read but **not** applied: a primer sitting in a
repeat is exactly the case worth reporting.

```bash
python ispcr_scan.py --genome .cache/genome/danRer11.2bit \
    --candidates primers_out/candidates.json -o primers_out/specificity.json
```

The `.2bit` (421 MB, `hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/`)
lives in `.cache/genome/` and is not in git.

**The parser is checked against an independent source.** A home-grown binary
reader that is subtly wrong -- endianness, the T/C/A/G code order, an off-by-one
in the N blocks -- yields sequence that still looks like DNA and is wrong
everywhere, which would quietly invalidate every specificity call. So four
regions on four chromosomes are compared byte for byte against UCSC's REST API
(different code path, different server, 1-based inclusive against the parser's
0-based). The scanner's own matching logic is unit-tested on synthetic sequence:
a planted amplicon is found at the right size, a 5' mismatch is tolerated, a
single 3'-seed mismatch is rejected, and a duplicated amplicon is reported
twice.

### Feeding it back into the design

`primer_design.py --specificity primers_out/specificity.json` ranks the
shortlist **specificity first**: a pair that amplifies a second locus does not
genotype anything, however good its Tm. A pair whose specificity was never
measured sorts between "clean" and "has an off-target" -- not promoted over a
pair known to be clean, not condemned on an absence of evidence. The result
travels into `primers.csv` (`ispcr_products`, `ispcr_off_target`, `ispcr_hits`),
the report's **In-silico PCR** column, and the conditions panel, which records
the genome, seed length and mismatch ceiling the numbers came from.

Keyed on the oligo pair, never on a rank: ranks are exactly what re-ranking
changes.

### Result on the test set

**332/332 candidate pairs across the twelve loci predict exactly one product in
the primary assembly**, the twenty-four ordered pairs included. Every ordered
pair's *intended* product was found at the expected chromosome and size, which
is the scan's own self-check: a pair whose intended product is missing means the
scanner is wrong, not the primer.

Specificity therefore did not reorder anything here -- it confirmed a choice
made on other grounds. That is a result, not a non-event: it is the difference
between twelve pairs that are probably fine and twelve pairs that have been
checked.

**The alt-contig trap, which nearly inverted this.** The first run of the scan
reported an off-target for **all twelve ordered yap1 pairs** -- and for 101 of
the 189 candidates then in the shortlist (172 of 332 in the final, deeper one). Every one of those hits was on `chr18_KZ115577v1_alt`, an
alternate-haplotype contig covering the yap1 region. danRer11 contains 26
primary sequences, 967 unplaced scaffolds and **930 alt contigs**; an alt contig
is a region of the primary assembly written down a second time, and a fish
carries one version of it, so a product predicted there is the intended product
seen twice, not a second band.

Had that been taken at face value, specificity-first ranking would have demoted
every yap1 primer for a reason that does not exist at the bench and replaced
them with worse ones. The scan now classifies each sequence (`seq_class`) and
counts only primary and unplaced hits, reporting alt-haplotype copies in a
separate field. Unplaced scaffolds still count: they are real sequence that
simply has no coordinate yet.

The lesson generalises past this gene. **A specificity check that cannot tell an
assembly duplicate from a second locus is worse than no check**, because it is
confidently wrong in the direction of changing a good answer.

## refresh_job.py -- re-point a report at a live CHOPCHOP job

CHOPCHOP expires finished jobs, after which every `details/<rank>` link in a
generated report 404s. The report's data is unaffected (it was captured into the
HTML at generation time); only the outbound links die.

```bash
python refresh_job.py --check out_geneA/geneA_danRer11.html   # dry run
python refresh_job.py --apply --job <id> out_geneA/geneA_danRer11.html
```

It submits a new job with the identical option array and matches old guides to
new ones by `(chrom, start, strand, protospacer)` -- **not** by rank, since rank
is the thing under suspicion. Only if every guide matches and every match kept
its rank does it rewrite anything, and then it rewrites exactly two strings,
`job` and `jobUrl`. Everything else is left alone because nothing is re-derived.

Why the rank guard matters: starred guides live in `localStorage` under
`gr-pins-<gene>-<build>`, which does not contain the job id, so they survive a
refresh -- but they are stored as *rank numbers*. Silently renumbering would
move a saved selection onto different guides. The script aborts instead.

Verified 2026-11-13: both genes matched 240/240 and 241/241 with zero
renumbering, so CHOPCHOP's ranking is reproducible for a fixed option array.
Diffing the payloads before and after confirmed `job` and `jobUrl` were the only
keys that changed, with the `guides` array byte-identical.

## The report's Primers tab

`primer_design.py` writes `candidates.json` alongside `primers.csv`
(`--max-candidates`, default 8). The report injects it and renders a second tab.

**Conditions panel.** The top of the tab states what the search actually ran
under: the shared Tm and how it was found, product range, minimum cut-to-end
distance, dTm ceiling, primer length and GC, which variant sets were masked,
fetch window, salt, Primer3 version and run time, plus the relaxation ladder in
order. `candidates.json` carries a `params` block for this -- without it the
report showed numbers with no way to tell what produced them.

Note `n_loci` in that block. The Tm scan optimises across **every locus in the
run**, which is usually more genes than one report shows, so the fraction has to
be reported against the recorded total. Using the page's own locus count printed
"best 11/6 loci".

**Re-run builder.** A collapsible panel with the run's settings pre-filled.
Changing anything assembles the exact `primer_design.py` command, emitting only
the flags that differ from the current run. The page is a static file and cannot
execute it -- a button that silently did nothing would be worse than no button
-- so it produces a command to copy and run in `crispr_tools/`.

**Ranking and selection.** The `#` column is the rank: pair 1 is what the
designer would order, sorted by dTm then Primer3 penalty. There is no separate
"chosen" marker, because the rank already says it. Each row has a checkbox and
**several pairs per guide can be ticked** -- ordering a backup pair is a second
checkbox. Everything ticked is drawn on that locus's track, each pair in its own
lane against the same guide and coordinates, and goes into the exports.

Selections persist per gene in `localStorage` under
`gr-primerpick-<gene>-<build>`, stored as `{sig, sel}`. `sig` is a per-locus
fingerprint of the shortlist and its recommendation, shipped in
`candidates.json`; when it no longer matches, the stored picks are dropped.
Without that, a pick made against an *earlier* design survived as a bare rank
number and silently masked the new recommendation -- and rank 2 after a
re-design is a different pair of oligos. Indices that no longer exist are
dropped too, and if that empties a locus it falls back to rank 1, so a locus is
never left silently unselected.

**Two pairs are pre-selected per guide, always.** `recommend_pairs()` in
`primer_design.py` picks rank 1 as the primary, then prefers a backup that binds
**neither** of the primary's sites. Primer3 usually holds one side fixed and
walks the other, so consecutive ranks often re-use a site -- a backup sharing
the reverse primer fails for the same reason the primary did.

**Sites, not sequences.** Primer3 routinely returns a primer that is its
neighbour extended or trimmed by a base: a different string, the same binding
site, the same SNP under it, the same mispriming partner. So independence is
decided on coordinates -- overlap of more than half the shorter oligo is the
same site -- not on string equality. The report applies the identical rule, so
the page and the designer cannot disagree about the word (they did, briefly:
the page compared `p.fwd === q.fwd` and so stayed silent on exactly the backups
that needed the warning).

Where no independent pair exists the near-duplicate is **still recommended and
still ordered**, flagged `rec_independent: false`, with a `rec_why` saying what
is and is not insured. Twelve pairs per gene is the plate layout, and a backup
that shares a failure mode is still a second shot at an oligo that was simply
synthesised badly. The locus track shows the warning per side and per kind
("share the forward primer and sit on the same reverse site").

That search needs depth. The first independent pair is routinely rank 11, 13,
19, 22 or 27, so Primer3 is asked for **40** pairs per locus (`--num-return`,
default 40) and `--max-candidates 40` keeps them: measured on these twelve loci,
40 gives **8/12** independent backups where 20 gives 7/12 and 60 gives 8/12
again. The deeper list costs one extra Primer3 call per locus and nothing else. The table still shows the top
8, with a "Show all N pairs" toggle, and any selected pair deeper than that
stays visible regardless -- a table must not claim a selection it does not show.

Result on the twelve-guide test set: **8/12 get a genuine independent
backup**; the other four are ordered with the flagged near-duplicate. The
hardest locus is the one where Primer3 returns only 3 pairs at any tier, so no
depth setting can help it.

`rec_independent` is a boolean in `candidates.json`, deliberately not left to be
inferred from `rec_why`: the prose for the negative case contains the word
"independent" ("no independent pair exists"), and string-matching on it reads
the exact opposite of the truth. That mistake was made once already.

Exports, all client-side from the current selection:

| button | what you get |
|---|---|
| Export selected pairs (CSV) | every ticked pair, 37 columns |
| Export all candidates (CSV) | every shortlisted pair, with an `is_selected` flag |
| Order sheet -- locus primers | Name/Sequence/Scale/Purification |
| Order sheet -- with Illumina overhangs | same, with the 33/34 nt tails prepended |
| Reset to rank 1 | clears the stored selection |

Order-sheet names gain a `_p<rank>` suffix only when more than one pair is
ordered for that locus, so the common case keeps a short tube label.

The CSV carries a `specificity` column spelling out that it was not checked. A
spreadsheet gets mailed around and separated from the page that warned about it.

Verified in-browser on both genes: 12 condition chips, the command builder emits
only changed flags, ticking pairs 3 and 5 gives `[1,3,5]` with three lanes drawn
and 8 selected CSV rows across the 6 loci, order names disambiguate only where
needed, selections survive a reload, and there are no console errors.
