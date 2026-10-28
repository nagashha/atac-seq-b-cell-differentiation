# Chromatin accessibility across B-cell differentiation

An ATAC-seq analysis comparing **bulk B cells with plasmablasts**, from raw reads through to transcription-factor motif enrichment. Built as a three-week internship project at the Department of Medical Bioinformatics, Universitätsmedizin Göttingen.

![Shell](https://img.shields.io/badge/shell-Bash-4EAA25)
![R](https://img.shields.io/badge/R-4.5.3-276DC3)
![DESeq2](https://img.shields.io/badge/DESeq2-1.50.2-8A2BE2)
![License](https://img.shields.io/badge/license-MIT-green)

**[Read the full report (PDF)](report/report.pdf)** · **[Discussion slides](slides/discussion_deck.pptx)** · Data: [GSE118189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118189) (Calderon et al. 2019, *Nature Genetics*)

---

## The question

Resting B cells and plasmablasts are two stages of one lineage. Getting from the first to the second means dismantling the B-cell identity programme and installing an antibody-secreting one. Which parts of the genome open and close as that happens, and which transcription factors drive it?

## Headline result

| | |
|---|---|
| Consensus regions quantified | 155,948 across 7 samples |
| Regions tested after filtering | 101,592 |
| Differentially accessible (FDR < 0.05, \|log2FC\| ≥ 1) | **14,662** |
| More accessible in plasmablasts | 7,976 |
| More accessible in bulk B | 6,686 |

<p align="center">
  <img src="report/figures/fig3c_pca.png" width="47%" alt="PCA of variance-stabilised accessibility, separating the two cell types along PC1">
  <img src="report/figures/fig4a_volcano.png" width="47%" alt="Volcano plot of differential accessibility">
</p>

PC1 accounts for 87% of the variance among the 500 most variable regions and separates the two cell types cleanly, with no donor structure visible.

## What the motifs say

Motif enrichment splits cleanly by direction, which is the central result of the project:

| More accessible in plasmablasts | *p* | | More accessible in bulk B | *p* |
|---|---|---|---|---|
| IRF4 | 1e-202 | | Elf4 (ETS) | 1e-324 |
| IRF3 | 1e-147 | | PU.1 | 1e-306 |
| IRF8 | 1e-136 | | ELF5 (ETS) | 1e-289 |
| bZIP:IRF (BATF) | 1e-111 | | SpiB | 1e-260 |

The results are consistent with an **ETS-associated programme giving way to an IRF-associated one** — PU.1 and SpiB being the factors that hold B-cell identity in place. Motif enrichment measures sequence composition, not occupancy, so this is a regulatory shift the data are compatible with rather than one they demonstrate.

Two analyses converge on it from different directions. IRF4 was picked as a target gene from promoter-proximal differential regions, and separately comes out as the top enriched motif. PU.1 appears both in HOMER's sequence scan and in g:Profiler's curated TF-target annotation — different methods, different databases, same factor. Both start from the same differential set, so they are convergent rather than independent.

## Pipeline

| Stage | Tools | Output |
|---|---|---|
| Manual pre-processing | Cutadapt, Bowtie2, SAMtools, Picard, BEDTools, MACS3 | 1,035 peaks; validated against an nf-core reference run |
| Consensus peak set | BEDTools, featureCounts | 155,948 regions × 7 samples |
| Differential accessibility | DESeq2, apeglm | 14,662 differential regions |
| Downstream | HOMER, g:Profiler, IGV | motif enrichment, GO terms, six target genes |

Week 1 was run by hand rather than through nf-core, deliberately, to understand what each filter actually does. Comparing the two pipelines on the same library gave duplication within three percentage points (76.59% vs 73.76%) and the same usable read fraction at the end of the pipeline (20.21% vs 20.73%).

## Decisions worth explaining

Four choices in the analysis were not obvious, and the reasoning behind them is the part of this project I would most want a reader to see.

**HOMER background.** Motif enrichment was run against the full consensus peak set rather than random genomic sequence. That changes the question from "what is in open chromatin" to "what distinguishes regions that changed from regions open in both cell types" — which is the one worth asking.

**GO background.** Enrichment against the whole annotated genome would include genes the experiment had no way of reporting. The defensible universe is the 16,776 protein-coding genes (of 19,419) with a TSS within 5 kb of any consensus region — a 13.6% reduction, which bounds how much the original choice could have mattered. [`results/gprofiler_background.txt`](results/gprofiler_background.txt) is that list, `scripts/03_annotation_and_motifs.sh` derives it, and [`docs/go_background.md`](docs/go_background.md) reports it.

**Design constraint.** The assignment specified a paired design (`~ Donor + cell_type`), but only 2 of 5 donors contributed both cell types. The main analysis therefore uses `~ cell_type` on all seven samples, with the paired design run on the four matched samples as a sensitivity check. It recovered 3,692 of its 3,696 differential regions from the main analysis, with complete agreement on direction. Two matched donors cannot rule out residual confounding, and the report says so.

**Filter threshold.** The slides specify a low-count filter of ≥ 10 counts in ≥ 3 samples; the assignment PDF specifies 4. Both were run to the end — see `scripts/R/05_filter_sensitivity.R` and [`docs/filter_sensitivity.md`](docs/filter_sensitivity.md) — so the choice rests on measured concordance rather than on an argument.

**Limitations are documented, not hidden.** Nearest-TSS gene assignment is an assumption rather than a measurement. Consensus regions are merges of broadPeak calls, so their centres are not called summits. A one-base BED-to-SAF offset was found after the fact, quantified, and corrected in the script here. All of it is in the report's discussion.

## Layout

```
report/     LaTeX source, figures, compiled PDF
scripts/    shell pipeline (weeks 1-3)
  R/        DESeq2 analysis and sensitivity checks
results/    differential regions, gene lists, motif tables, GO background
slides/     discussion-session deck
docs/       follow-up analyses done after submission
```

Paths in the scripts point at a university cluster and are set through environment variables (`PEAK_DIR`, `BAM_DIR`, `HOMER_REF`). Raw FASTQ and BAM files are not included; they are available from GEO under GSE118189.

To rebuild the report:

```bash
cd report && latexmk -pdf report.tex
```

## Software

Cutadapt 4.7 · Bowtie2 2.5.5 · SAMtools 1.24 · Picard 3.4.0 · BEDTools 2.31.1 · MACS3 3.0.2 · MultiQC 1.35 · Subread/featureCounts 2.1.1 · R 4.5.3 · DESeq2 1.50.2 · apeglm 1.32.0 · HOMER 4.11 · g:Profiler e114_eg62_p19

## Note

An AI assistant was used for language refinement, figure preparation and technical troubleshooting. The analytical decisions, interpretation and conclusions are my own.
