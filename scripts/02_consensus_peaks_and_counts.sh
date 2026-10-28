#!/usr/bin/env bash
#
# Week 2, days 1-2 : build a consensus peak set across samples and quantify it.
#
# Peak calls do not line up between samples, so accessibility can only be
# compared after every sample has been measured over the same set of regions.
# This script pools the per-sample calls, merges overlaps into one consensus
# set, and counts fragments per region per sample.
#
# Inputs   nf-core/atacseq per-sample broadPeak calls and filtered BAMs
# Outputs  consensus_peaks.bed, consensus_peaks.saf, counts_clean.txt
# Requires bedtools, featureCounts (subread)
#
set -euo pipefail

# --------------------------------------------------------------- settings ---
PEAK_DIR="${PEAK_DIR:-./data/peaks}"      # *_peaks.broadPeak
BAM_DIR="${BAM_DIR:-./data/bams}"         # *.sorted.bam (+ .bai)
OUT_DIR="${OUT_DIR:-./week2}"
THREADS="${THREADS:-4}"

# Filename patterns, kept in one place so the script survives a different
# nf-core naming convention.
PEAK_SUFFIX="_REP1.mLb.clN_peaks.broadPeak"
BAM_SUFFIX="_REP1.mLb.clN.sorted.bam"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# Derive the sample list from the files themselves rather than typing it out,
# so the list cannot drift from what is actually present on disk.
mapfile -t SAMPLES < <(
  find "$PEAK_DIR" -maxdepth 1 -name "*${PEAK_SUFFIX}" -printf '%f\n' |
  sed "s/${PEAK_SUFFIX}\$//" | sort
)
(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: no peak files in $PEAK_DIR" >&2; exit 1; }
printf 'samples (%d): %s\n' "${#SAMPLES[@]}" "${SAMPLES[*]}"

# ------------------------------------------- 1. consensus peak set ----------
# Pool coordinates from every sample, sort, merge overlaps. Peaks called in
# only one sample are kept deliberately: the other samples then receive
# near-zero counts there, and that near-zero is itself the measurement.
for sample in "${SAMPLES[@]}"; do
    cut -f1-3 "${PEAK_DIR}/${sample}${PEAK_SUFFIX}"
done | sort -k1,1 -k2,2n | bedtools merge -i - > consensus_peaks.bed

printf 'consensus regions: %s\n' "$(wc -l < consensus_peaks.bed)"

# Width distribution. Worth checking, because merging overlapping calls can in
# principle chain into very wide regions, and region width feeds both the
# counts and the HOMER window used in script 03.
awk '{print $3-$2}' consensus_peaks.bed | sort -n | awk '
  {w[NR]=$1; s+=$1}
  END {
    printf "width  min %d  median %d  mean %.0f  p95 %d  p99 %d  max %d\n",
           w[1], w[int((NR+1)/2)], s/NR,
           w[int(0.95*NR)], w[int(0.99*NR)], w[NR]
  }'

# ------------------------------------------- 2. BED -> SAF ------------------
# BED is 0-based half-open; SAF is 1-based inclusive at both ends. The start
# therefore needs +1, or every counting interval is one base wider on the left.
#
# The submitted analysis omitted this +1. The offset is identical for every
# region and every sample, so it cannot plausibly have moved the between-group
# comparison, but it was a real error and the corrected form is what runs here.
awk 'BEGIN { OFS="\t"; print "GeneID","Chr","Start","End","Strand" }
     { print "peak_"NR, $1, $2+1, $3, "." }' \
     consensus_peaks.bed > consensus_peaks.saf

# ------------------------------------------- 3. fragment counting -----------
# Build the BAM list in the SAME order as SAMPLES. If these two ever diverge,
# the count columns carry the wrong sample labels and nothing downstream can
# detect it, because the metadata is built from the same (wrong) header.
BAMS=()
for sample in "${SAMPLES[@]}"; do
    bam="${BAM_DIR}/${sample}${BAM_SUFFIX}"
    [[ -f "$bam"       ]] || { echo "ERROR: missing BAM: $bam" >&2; exit 1; }
    [[ -f "${bam}.bai" ]] || { echo "ERROR: missing index: ${bam}.bai" >&2; exit 1; }
    BAMS+=("$bam")
done

# -p --countReadPairs  count fragments, not individual reads
# -B                   require both ends of a pair to be successfully aligned
#                      (this is NOT a requirement that both ends fall in the
#                       same feature, which is what the first draft claimed)
# -C                   discard chimeric pairs, i.e. ends on different chromosomes
featureCounts \
  -F SAF -a consensus_peaks.saf \
  -o counts_matrix.txt \
  -p --countReadPairs -B -C -T "$THREADS" \
  "${BAMS[@]}"

cat counts_matrix.txt.summary

# ------------------------------------------- 4. tidy the matrix -------------
# featureCounts layout: 1=GeneID 2=Chr 3=Start 4=End 5=Strand 6=Length 7+=samples
# The header is rebuilt from SAMPLES so the labels cannot drift from the
# columns; featureCounts would otherwise name them after full BAM paths.
header=$(IFS=$'\t'; echo "PeakID"$'\t'"${SAMPLES[*]}")

awk -v hdr="$header" '
     NR==1 && /^#/ { next }                       # featureCounts command line
     NR==2         { print hdr; next }            # replace its column names
     { out=$1; for (i=7; i<=NF; i++) out=out"\t"$i; print out }
' counts_matrix.txt > counts_clean.txt

# Sanity check: one header line plus one line per consensus region.
expected=$(( $(wc -l < consensus_peaks.bed) + 1 ))
actual=$(wc -l < counts_clean.txt)
[[ "$expected" -eq "$actual" ]] ||
  { echo "ERROR: counts_clean.txt has $actual lines, expected $expected" >&2; exit 1; }

echo "counts_clean.txt: $actual lines (1 header + $(( actual - 1 )) regions)"
head -3 counts_clean.txt
