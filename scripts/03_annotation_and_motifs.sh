#!/usr/bin/env bash
#
# Week 3 : peak-to-gene annotation, target genes, HOMER motif enrichment,
#          and the gene universe used as the GO enrichment background.
#
# Turns a list of differential coordinates into biology: which genes sit
# nearby, and which transcription-factor motifs distinguish the regions that
# changed from the regions that stayed open in both cell types.
#
# Inputs   DARs_significant.tsv and consensus_peaks.saf from week 2
# Outputs  DAR_{up,down}*.bed, gene lists, gprofiler_background.txt,
#          homer_{up,down}/ and the top-motif tables
# Requires bedtools, HOMER, wget
#
set -euo pipefail

# --------------------------------------------------------------- settings ---
WEEK2_DIR="${WEEK2_DIR:-./week2}"
OUT_DIR="${OUT_DIR:-./week3}"
THREADS="${THREADS:-4}"
GTF_URL="${GTF_URL:-https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz}"

# Promoter window for calling a peak-to-gene assignment credible, and the
# window used to define the GO background. Same number by design: the
# background is "genes this experiment could have reported".
PROMOTER_BP="${PROMOTER_BP:-5000}"

# HOMER's own genome packages were not installable in the shared environment
# (no write permission on the shared HOMER directory), so a GRCh38
# primary-assembly FASTA is passed to findMotifsGenome.pl directly.
HOMER_REF="${HOMER_REF:?set HOMER_REF to a GRCh38 primary-assembly FASTA}"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
cp "${WEEK2_DIR}/DARs_significant.tsv" "${WEEK2_DIR}/consensus_peaks.saf" .

# ------------------------------------------- 1. split DARs by direction -----
# DARs_significant.tsv columns:
#   1=baseMean 2=log2FoldChange 3=lfcSE 4=pvalue 5=padj 6=PeakID
# ($2+0) forces a numeric comparison. Without it awk compares the
# scientific-notation strings lexically and silently returns nothing.
split_dars () {               # $1 = comparison, $2 = output file
  awk -v op="$1" '
    NR==FNR && FNR>1 { coord[$1] = $2"\t"$3"\t"$4; next }
    FNR>1 && $6 in coord {
      lfc = $2 + 0
      if ((op == "up"   && lfc >=  1) ||
          (op == "down" && lfc <= -1)) print coord[$6]"\t"$6"\t0\t+"
    }' consensus_peaks.saf DARs_significant.tsv > "$2"
}
split_dars up   DAR_up_plasmablast.bed
split_dars down DAR_down_bulkB.bed

# HOMER background: every consensus region, i.e. the whole accessible universe.
awk 'FNR>1 { print $2"\t"$3"\t"$4"\t"$1"\t0\t+" }' \
    consensus_peaks.saf > consensus_all.bed

wc -l DAR_up_plasmablast.bed DAR_down_bulkB.bed consensus_all.bed

# ------------------------------------------- 2. motif enrichment ------------
# -bg consensus_all.bed  compare the regions that changed against all open
#                        regions rather than against random genomic sequence.
#                        This changes the question from "what is in open
#                        chromatin" to "what distinguishes what changed".
# -size 200              +/- 100 bp around the CENTRE of each merged consensus
#                        interval. This centre is not a called summit: the
#                        intervals are merges of broadPeak calls across
#                        samples. With a median consensus width of 350 bp the
#                        window covers most of a typical region, but for the
#                        widest regions it may miss the true signal maximum.
run_homer () {                # $1 = target bed, $2 = output dir, $3 = log
  findMotifsGenome.pl "$1" "$HOMER_REF" "$2" \
    -bg consensus_all.bed -size 200 -p "$THREADS" > "$3" 2>&1
}
run_homer DAR_up_plasmablast.bed homer_up/   homer_up.log   & PID_UP=$!
run_homer DAR_down_bulkB.bed     homer_down/ homer_down.log & PID_DOWN=$!

# ------------------------------------------- 3. gene annotation -------------
# Runs while HOMER works in the background.
GTF=$(basename "$GTF_URL")
[[ -f "$GTF" ]] || wget -q "$GTF_URL"

# The TSS depends on strand: gene start on '+', gene END on '-'. Getting this
# backwards would misplace every minus-strand TSS, which is half the genome,
# and would do so without any error.
zcat "$GTF" | awk -F'\t' '
  $3 == "gene" {
    tss = ($7 == "+") ? $4 : $5
    name = ""; type = ""
    if (match($9, /gene_name "[^"]+"/))    name = substr($9, RSTART+11, RLENGTH-12)
    if (match($9, /gene_biotype "[^"]+"/)) type = substr($9, RSTART+14, RLENGTH-15)
    if (name != "" && type == "protein_coding")
      print $1"\t"(tss-1)"\t"tss"\t"name"\t0\t"$7
  }' | sort -k1,1 -k2,2n > tss_protein_coding.bed

printf 'protein-coding TSS: %s\n' "$(wc -l < tss_protein_coding.bed)"

sort -k1,1 -k2,2n DAR_up_plasmablast.bed > DAR_up.sorted.bed
sort -k1,1 -k2,2n DAR_down_bulkB.bed     > DAR_down.sorted.bed

# Nearest protein-coding TSS to each DAR, with the distance in column 13.
# DAR columns 1-6, nearest TSS columns 7-12, distance column 13.
bedtools closest -a DAR_up.sorted.bed   -b tss_protein_coding.bed -d > DAR_up_annotated.txt
bedtools closest -a DAR_down.sorted.bed -b tss_protein_coding.bed -d > DAR_down_annotated.txt

# ------------------------------------------- 4. gene lists for GO -----------
# The distance cutoff matters: a peak 200 kb from a TSS may well loop to an
# entirely different gene, so only promoter-proximal assignments go forward.
for dir in up down; do
  awk -v d="$PROMOTER_BP" '$13 <= d { print $10 }' "DAR_${dir}_annotated.txt" |
    grep -vE '^\.$|^$' | sort -u > "genes_${dir}_promoter.txt"
done
wc -l genes_up_promoter.txt genes_down_promoter.txt

# ------------------------------------------- 5. GO background ---------------
# Enrichment against the whole annotated genome includes genes this experiment
# had no way of reporting. The defensible universe is every protein-coding
# gene whose TSS lies within the same promoter window of ANY consensus region.
awk 'FNR>1 { print $2"\t"$3-1"\t"$4 }' consensus_peaks.saf |
  sort -k1,1 -k2,2n > consensus_regions.bed

bedtools window -a tss_protein_coding.bed -b consensus_regions.bed \
                -w "$PROMOTER_BP" -u |
  cut -f4 | sort -u > gprofiler_background.txt

printf 'GO background: %s of %s protein-coding genes\n' \
  "$(wc -l < gprofiler_background.txt)" "$(cut -f4 tss_protein_coding.bed | sort -u | wc -l)"

# Confirm the query lists are subsets of the background, which they must be.
for dir in up down; do
  stray=$(comm -23 "genes_${dir}_promoter.txt" gprofiler_background.txt | wc -l)
  [[ "$stray" -eq 0 ]] ||
    echo "WARNING: $stray genes in genes_${dir}_promoter.txt are outside the background" >&2
done

# ------------------------------------------- 6. marker sanity check ---------
# Check which canonical markers actually turn up, BEFORE choosing targets.
echo "--- plasma-cell programme (expected in the up list) ---"
grep -wE "PRDM1|XBP1|IRF4|CD38|SDC1|MZB1|TNFRSF17" genes_up_promoter.txt || true
echo "--- B-cell identity (expected in the down list) ---"
grep -wE "PAX5|BACH2|MS4A1|CD19|CR2|CIITA|EBF1"    genes_down_promoter.txt || true

# ------------------------------------------- 7. collect HOMER --------------
# Wait before reading the result files, or head runs on a half-written table.
wait "$PID_UP" "$PID_DOWN"
for dir in up down; do
  head -16 "homer_${dir}/knownResults.txt" | cut -f1,3,5,7,9 > "homer_${dir}_top15.tsv"
done

# GO enrichment itself was run through the g:Profiler web interface
# (https://biit.cs.ut.ee/gprofiler/gost), organism Homo sapiens, with
# gprofiler_background.txt pasted into "Custom background" and the three
# query lists being genes_up_promoter.txt, genes_down_promoter.txt and the
# six target genes.
