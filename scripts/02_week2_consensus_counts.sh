#!/usr/bin/env bash
#
# Week 2, days 1-2 - build a consensus peak set across samples and quantify it.
#
# Individual peak calls do not line up between samples, so accessibility can
# only be compared after every sample is measured over the same regions.
#
# Requires: bedtools, samtools, featureCounts (subread)
#
set -euo pipefail

# ---------------------------------------------------------------- settings --
# Point these at your own copies of the nf-core/atacseq output.
PEAK_DIR="${PEAK_DIR:-./data/peaks}"   # *_peaks.broadPeak files
BAM_DIR="${BAM_DIR:-./data/bams}"      # *.sorted.bam files (+ .bai)
OUT_DIR="${OUT_DIR:-./week2}"
THREADS=4

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# Build the sample list from the files themselves rather than typing it out.
SAMPLES=($(ls "$PEAK_DIR"/*_peaks.broadPeak | xargs -n1 basename | sed 's/_REP1.*//'))
echo "Samples (${#SAMPLES[@]}): ${SAMPLES[*]}"

# ------------------------------------------------ 1. consensus peak set -----
# Pool coordinates from every sample, sort, merge overlaps. Peaks unique to a
# single sample are kept: other samples simply receive near-zero counts there,
# and that near-zero is itself informative.
for sample in "${SAMPLES[@]}"; do
    awk '{print $1"\t"$2"\t"$3}' "${PEAK_DIR}/${sample}_REP1.mLb.clN_peaks.broadPeak"
done | sort -k1,1 -k2,2n | bedtools merge -i - > consensus_peaks.bed

wc -l consensus_peaks.bed

# ------------------------------------------------ 2. BED -> SAF -------------
# NOTE: BED is 0-based half-open, SAF is 1-based inclusive, so the start
# should be $2+1. The line below reproduces the conversion as it was actually
# run, which expands every interval by one base on the left. The offset is
# systematic across all peaks and samples, so it does not bias the
# between-group comparison, but the corrected form is:
#     {print "peak_"NR,$1,$2+1,$3,"."}
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"}
     {print "peak_"NR,$1,$2,$3,"."}' consensus_peaks.bed > consensus_peaks.saf

# ------------------------------------------------ 3. fragment counting ------
# Build the BAM list in the SAME order as SAMPLES, or the count columns will
# be mislabelled downstream.
BAMS=()
for sample in "${SAMPLES[@]}"; do
    BAMS+=("${BAM_DIR}/${sample}_REP1.mLb.clN.sorted.bam")
done

for b in "${BAMS[@]}"; do
    [[ -f "$b"       ]] || { echo "ERROR: BAM not found: $b";   exit 1; }
    [[ -f "${b}.bai" ]] || { echo "ERROR: index not found: $b"; exit 1; }
done

# -p --countReadPairs  count fragments, not individual reads
# -B                   count only fragments with both ends aligned
# -C                   discard chimeric pairs (ends on different chromosomes)
featureCounts \
  -F SAF -a consensus_peaks.saf \
  -o counts_matrix.txt \
  -p --countReadPairs -B -C -T "$THREADS" \
  "${BAMS[@]}"

cat counts_matrix.txt.summary

# ------------------------------------------------ 4. tidy the matrix --------
# featureCounts columns: 1=GeneID 2=Chr 3=Start 4=End 5=Strand 6=Length 7+=samples
# Header is rebuilt from SAMPLES so the labels cannot drift from the columns.
HEADER="PeakID"
for s in "${SAMPLES[@]}"; do HEADER="${HEADER}\t${s}"; done

awk -v hdr="$HEADER" '
     NR==1 && /^#/ {next}
     NR==2 { print hdr; next }
     { line=$1; for (i=7; i<=NF; i++) line = line "\t" $i; print line }' \
     counts_matrix.txt > counts_clean.txt

head -3 counts_clean.txt
wc -l counts_clean.txt
