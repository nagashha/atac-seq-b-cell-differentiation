#!/usr/bin/env bash
#
# Week 1 - manual ATAC-seq pre-processing of a single paired-end library.
#
# Runs every step by hand rather than through nf-core/atacseq, so that the
# purpose of each filter is explicit. Input is a chromosome-22 subset of
# SRR7650763 (an immature NK-cell library, used here purely as a methods
# exercise - see the report, section 1.2).
#
# Requires: cutadapt, bowtie2, samtools, picard, bedtools, macs3, multiqc
#
set -euo pipefail

# ---------------------------------------------------------------- settings --
SAMPLE=SRR7650763
REF_INDEX=ref/GRCh38_chr22          # bowtie2 index prefix (no extension)
BLACKLIST=hg38-blacklist.ensembl.bed
THREADS=2                           # fixed at 2 for the hardware benchmark
EFFECTIVE_GENOME=4.0e7              # chr22 mappable size, NOT whole-genome 'hs'

mkdir -p output_cutadapt macs2_out

# ------------------------------------------------ 1. adapter trimming -------
# ATAC-seq libraries carry the Tn5 mosaic-end adapter; tagmentation is
# symmetric, so the same sequence applies to both mates.
cutadapt \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -o output_cutadapt/trimmed_R1.fastq.gz \
  -p output_cutadapt/trimmed_R2.fastq.gz \
  "${SAMPLE}_1_chr22.raw.fastq.gz" \
  "${SAMPLE}_2_chr22.raw.fastq.gz" \
  > output_cutadapt/cutadapt.log

# ------------------------------------------------ 2. genome alignment -------
# -X 2000 is the one ATAC-specific setting. Bowtie2's default 500 bp insert
# limit would flag every fragment above it as improperly paired, and the
# -f 2 filter below would then discard them: 8.8% of this library.
bowtie2 -p "$THREADS" -X 2000 \
  -x "$REF_INDEX" \
  -1 output_cutadapt/trimmed_R1.fastq.gz \
  -2 output_cutadapt/trimmed_R2.fastq.gz \
  -S sample.sam \
  2> "${SAMPLE}_bowtie2.log"

# ------------------------------------- 3. filtering and duplicate removal ---
# -q 30   MAPQ >= 30: drops low-confidence and multi-mapping alignments
# -f 2    keep properly paired reads only
# -F 1804 drop unmapped(4) + mate unmapped(8) + secondary(256)
#         + vendor QC fail(512) + flagged duplicate(1024)
samtools view -b -q 30 -F 1804 -f 2 sample.sam > sample.filtered.bam
samtools sort -@ "$THREADS" -o sample.sorted.bam sample.filtered.bam

# Picard requires an @RG tag on every alignment; bowtie2 does not add one,
# and MarkDuplicates fails with a NullPointerException without it.
picard AddOrReplaceReadGroups \
  I=sample.sorted.bam O=sample.rg.bam \
  RGID=1 RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 "RGSM=${SAMPLE}"

picard MarkDuplicates \
  I=sample.rg.bam O=sample.dedup.bam \
  "M=${SAMPLE}.dup_metrics.txt" \
  REMOVE_DUPLICATES=true

samtools index sample.dedup.bam

# ------------------------------------------------ 4. blacklist filtering ----
# Verify the genome build first: a coordinate-based filter applied across
# builds fails silently and masks the wrong regions.
samtools view -H sample.dedup.bam | grep '^@SQ'   # expect SN:22 LN:50818468
head -3 "$BLACKLIST"                              # expect Ensembl naming

bedtools intersect -v -abam sample.dedup.bam -b "$BLACKLIST" > sample.clean.bam
samtools index sample.clean.bam

# ------------------------------------------------ 5. peak calling -----------
# -f BAMPE      use real fragment coordinates from mate pairs
# -g            effective genome size: chr22 only, not whole genome
# --keep-dup all duplicates already removed by Picard; do not filter twice
macs3 callpeak \
  -t sample.clean.bam \
  -f BAMPE \
  -g "$EFFECTIVE_GENOME" \
  -n "$SAMPLE" \
  --broad --broad-cutoff 0.1 \
  --keep-dup all \
  --outdir macs2_out \
  2> macs3.log

# ------------------------------------------------ 6. QC ---------------------
# Insert-size distribution is the single most informative ATAC-seq QC metric.
picard CollectInsertSizeMetrics \
  I=sample.dedup.bam \
  "O=${SAMPLE}.insert_size_metrics.txt" \
  "H=${SAMPLE}.insert_size_histogram.pdf" \
  M=0.5

multiqc . -o multiqc_out --force

# ------------------------------------------------ 7. summary numbers --------
echo "reads after dedup:     $(samtools view -c sample.dedup.bam)"
echo "reads after blacklist: $(samtools view -c sample.clean.bam)"
echo "peaks called:          $(wc -l < "macs2_out/${SAMPLE}_peaks.broadPeak")"

# FRiP - fraction of reads in peaks
IN=$(bedtools intersect -a sample.clean.bam \
       -b "macs2_out/${SAMPLE}_peaks.broadPeak" -u | samtools view -c)
TOT=$(samtools view -c sample.clean.bam)
echo "FRiP:                  $(echo "scale=4; $IN/$TOT" | bc)"
