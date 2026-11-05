# chopchop_pipeline

End-to-end CRISPR guide selection for zebrafish F0 knockouts, following
Kroll & Rihel's protocols. Submits (or reuses) a CHOPCHOP job, pulls everything
CHOPCHOP computes, annotates it against the RefSeq gene model, checks strain SNPs,
and writes an annotated CSV plus an interactive HTML report.

## What is in this folder

```
chopchop_pipeline.py    the script
report_template.html    REQUIRED - must sit next to the script
README.md               this file
```

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
- click any row for primers, full inDelphi stats, position and SNP detail
- **CC** opens that guide on CHOPCHOP; **SNP** re-queries SNPfisher live in a new tab
  (a real form POST, so it works for rows the script did not pre-query); **⧉** copies
  the 20 nt without the PAM, ready for IDT
- a **wild-type strain** selector (any / TL / WIK / fli1a-eGFP) that switches the
  SNP column, filter, scoring penalty and summary to that line's allele
  frequencies. `0 /2` means zero SNPs that matter in your strain, two in total
- export the filtered set back out to CSV

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

This is a convenience for exploring trade-offs, not a published metric. The underlying
CHOPCHOP rank is kept in its own column.

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
