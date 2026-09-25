# Chromatin accessibility across B-cell differentiation

An ATAC-seq analysis comparing **bulk B cells with plasmablasts**, from raw
reads through to transcription-factor motif enrichment. Built as a three-week
internship project at the Department of Medical Bioinformatics,
Universitätsmedizin Göttingen.

Data: [Calderon et al. 2019, *Nature Genetics*](https://doi.org/10.1038/s41588-019-0505-9) (GEO: GSE118189)

📄 **[Read the full report (PDF)](report/report.pdf)** · 📊 **[Slides](slides/discussion_deck.pptx)**

---

## The question

Resting B cells and plasmablasts are two stages of one lineage. Getting from
the first to the second means dismantling the B-cell identity programme and
installing an antibody-secreting one. Which parts of the genome open and close
as that happens, and which transcription factors drive it?

## Main results

**14,662 of 101,592 tested regions are differentially accessible** (FDR < 0.05,
|log2FC| ≥ 1) — 7,976 opening in plasmablasts and 6,686 closing. PC1 of the
500 most variable peaks accounts for 87% of the variance and separates the two
cell types cleanly.

Motif enrichment shows a clean split in the regions that changed:

| Opening in plasmablasts | p | | Staying open in bulk B | p |
|---|---|---|---|---|
| IRF4 | 1e-202 | | Elf4 (ETS) | 1e-324 |
| IRF3 | 1e-147 | | **PU.1** | 1e-306 |
| IRF8 | 1e-136 | | ELF5 (ETS) | 1e-289 |
| bZIP:IRF (BATF) | 1e-111 | | **SpiB** | 1e-260 |

The picture is an ETS-associated regulatory programme (PU.1, SpiB — the factors
that hold B-cell identity in place) giving way to an IRF-associated one.

Three analyses converge without being designed to. IRF4 was selected as a
target gene from promoter-proximal DARs, and independently turns out to be the
top enriched motif. PU.1 was found by HOMER scanning DNA sequence and again by
g:Profiler through curated TF–target annotation — different methods, different
databases, same factor.

<p align="center">
  <img src="report/figures/fig4a_volcano.png" width="45%">
  <img src="report/figures/fig3c_pca.png" width="45%">
</p>

## What I did

| Stage | Tools | What it produced |
|---|---|---|
| Manual pre-processing | Cutadapt, Bowtie2, SAMtools, Picard, BEDTools, MACS3 | 1,035 peaks; validated against an nf-core reference run |
| Consensus peak set | BEDTools, featureCounts | 155,948 regions × 7 samples |
| Differential accessibility | DESeq2, apeglm | 14,662 DARs |
| Downstream | HOMER, g:Profiler, IGV | motif enrichment, GO terms, six target genes |

Week 1 was deliberately run by hand rather than through nf-core, to understand
what each filter does. Comparing the two pipelines on the same library gave
duplication within three percentage points (76.59% vs 73.76%) and the same
usable read fraction at the end of the pipeline (20.21% vs 20.73%).

## Things worth noting

A few decisions in the report that were not obvious:

**HOMER background.** Motif enrichment was run against the full consensus peak
set rather than random genomic sequence. This changes the question from "what
is in open chromatin" to "what distinguishes regions that *changed* from
regions open in both cell types" — which is the biologically meaningful one.

**Design constraint.** The assignment specified a paired design
(`~ Donor + cell_type`), but only 2 of 5 donors contributed both cell types.
The main analysis therefore uses `~ cell_type` on all seven samples, with the
paired design run on the four matched samples as a sensitivity check. It
recovered 3,692 of its 3,696 DARs from the main analysis with complete
agreement on direction.

**Fragment-size distribution.** The insert-size histogram is multi-modal
(~45 bp nucleosome-free, ~210 bp mono-nucleosomal, ~354 bp di-nucleosomal),
a direct physical readout of chromatin packaging.

**Limitations are documented, not hidden.** Nearest-TSS gene assignment is an
assumption rather than a measurement; GO enrichment used a default rather than
a peak-derived background; a one-base BED-to-SAF offset was found after the
fact and its (systematic, non-biasing) effect quantified. All are set out in
the report's discussion.

## Repository layout

```
report/          LaTeX source, figures, compiled PDF
scripts/         shell pipeline (weeks 1–3)
  R/             DESeq2 analysis
results/         DAR tables, gene lists, motif results
slides/          discussion-session deck
```

Paths in the scripts point at a university cluster and are set through
environment variables (`PEAK_DIR`, `BAM_DIR`, `HOMER_REF`). Raw FASTQ and BAM
files are not included — they are available from
[GEO GSE118189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118189).

To rebuild the report:

```bash
cd report && latexmk -pdf report.tex
```

## Software

Cutadapt 4.7 · Bowtie2 2.5.5 · SAMtools 1.24 · Picard 3.4.0 · BEDTools 2.31.1 ·
MACS3 3.0.2 · MultiQC 1.35 · Subread/featureCounts 2.1.1 · R 4.5.3 ·
DESeq2 1.50.2 · apeglm 1.32.0 · HOMER 4.11 · g:Profiler (e114_eg62_p19)

## Note

An AI assistant was used for language refinement, figure preparation and
technical troubleshooting. The analytical decisions, interpretation and
conclusions are my own.
