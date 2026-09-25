#!/usr/bin/env bash
#
# Week 3 - peak-to-gene annotation, target genes, HOMER motif enrichment.
#
# Turns 14,662 coordinates into biology: which genes are nearby, which
# transcription factors bind the regions that changed.
#
# Requires: bedtools, samtools, HOMER
#
set -euo pipefail

# ---------------------------------------------------------------- settings --
WEEK2_DIR="${WEEK2_DIR:-./week2}"
# HOMER genome packages were not installable in the shared environment, so the
# reference FASTA is passed to findMotifsGenome.pl directly.
HOMER_REF="${HOMER_REF:?set HOMER_REF to a GRCh38 primary-assembly FASTA}"
THREADS=4

mkdir -p week3 && cd week3
cp "${WEEK2_DIR}/DARs_significant.tsv" "${WEEK2_DIR}/consensus_peaks.saf" .

# ------------------------------------------------ 1. split DARs by direction
# DARs_significant.tsv columns: 1=baseMean 2=log2FoldChange 3=lfcSE 4=pvalue
#                               5=padj     6=PeakID
# ($2+0) forces numeric comparison: scientific notation compares as text otherwise.
awk 'NR==FNR && FNR>1 {c[$1]=$2"\t"$3"\t"$4; next}
     FNR>1 && ($2+0)>=1  && $6 in c {print c[$6]"\t"$6"\t0\t+"}' \
     consensus_peaks.saf DARs_significant.tsv > DAR_up_plasmablast.bed

awk 'NR==FNR && FNR>1 {c[$1]=$2"\t"$3"\t"$4; next}
     FNR>1 && ($2+0)<=-1 && $6 in c {print c[$6]"\t"$6"\t0\t+"}' \
     consensus_peaks.saf DARs_significant.tsv > DAR_down_bulkB.bed

# Background for HOMER: the full consensus set
awk 'FNR>1 {print $2"\t"$3"\t"$4"\t"$1"\t0\t+"}' consensus_peaks.saf > consensus_all.bed

wc -l DAR_up_plasmablast.bed DAR_down_bulkB.bed consensus_all.bed

# ------------------------------------------------ 2. motif enrichment -------
# -bg consensus_all.bed : compare changed regions against ALL open regions,
#                         not random genomic background. This changes the
#                         question from "what is in open chromatin" to "what
#                         distinguishes regions that changed".
# -size 200             : +/-100 bp around the midpoint of each merged
#                         consensus interval. Note this midpoint is not a
#                         called summit, since intervals are merges of
#                         broadPeak calls across samples.
findMotifsGenome.pl DAR_up_plasmablast.bed "$HOMER_REF" homer_up/ \
  -bg consensus_all.bed -size 200 -p "$THREADS" > homer_up.log 2>&1 &
PID_UP=$!

findMotifsGenome.pl DAR_down_bulkB.bed "$HOMER_REF" homer_down/ \
  -bg consensus_all.bed -size 200 -p "$THREADS" > homer_down.log 2>&1 &
PID_DOWN=$!

# ------------------------------------------------ 3. gene annotation --------
# Runs while HOMER works in the background.
wget -q https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

# TSS depends on strand: gene start on '+', gene END on '-'. Getting this
# backwards would misplace every minus-strand TSS, i.e. half the genome.
zcat Homo_sapiens.GRCh38.110.gtf.gz | awk -F'\t' '
  $3=="gene" {
    tss = ($7=="+") ? $4 : $5
    name=""; type=""
    if (match($9, /gene_name "[^"]+"/))    { name=substr($9,RSTART+11,RLENGTH-12) }
    if (match($9, /gene_biotype "[^"]+"/)) { type=substr($9,RSTART+14,RLENGTH-15) }
    if (name!="" && type=="protein_coding")
      print $1"\t"(tss-1)"\t"tss"\t"name"\t0\t"$7
  }' | sort -k1,1 -k2,2n > tss_protein_coding.bed

sort -k1,1 -k2,2n DAR_up_plasmablast.bed > DAR_up.sorted.bed
sort -k1,1 -k2,2n DAR_down_bulkB.bed     > DAR_down.sorted.bed

bedtools closest -a DAR_up.sorted.bed   -b tss_protein_coding.bed -d > DAR_up_annotated.txt
bedtools closest -a DAR_down.sorted.bed -b tss_protein_coding.bed -d > DAR_down_annotated.txt
# output: DAR cols 1-6, nearest TSS cols 7-12, distance col 13

# ------------------------------------------------ 4. gene lists for GO ------
# 5 kb cutoff: a peak 200 kb away may loop to an entirely different gene.
awk '$13<=5000 {print $10}' DAR_up_annotated.txt   | sort -u > genes_up_promoter.txt
awk '$13<=5000 {print $10}' DAR_down_annotated.txt | sort -u > genes_down_promoter.txt

for f in genes_*.txt; do grep -vE '^\.$|^$' "$f" > tmp && mv tmp "$f"; done
wc -l genes_up_promoter.txt genes_down_promoter.txt

# Check which canonical markers actually turn up, before choosing targets.
echo "--- plasma-cell programme (expected in the up list) ---"
grep -wE "PRDM1|XBP1|IRF4|CD38|SDC1|MZB1|TNFRSF17" genes_up_promoter.txt || true
echo "--- B-cell identity (expected in the down list) ---"
grep -wE "PAX5|BACH2|MS4A1|CD19|CR2|CIITA|EBF1"    genes_down_promoter.txt || true

# ------------------------------------------------ 5. collect HOMER ----------
wait "$PID_UP" "$PID_DOWN"
head -16 homer_up/knownResults.txt   | cut -f1,3,5,7,9 > homer_up_top15.tsv
head -16 homer_down/knownResults.txt | cut -f1,3,5,7,9 > homer_down_top15.tsv

# GO enrichment was run through the g:Profiler web interface
# (https://biit.cs.ut.ee/gprofiler/gost) on genes_up_promoter.txt,
# genes_down_promoter.txt and the six target genes, organism Homo sapiens.
