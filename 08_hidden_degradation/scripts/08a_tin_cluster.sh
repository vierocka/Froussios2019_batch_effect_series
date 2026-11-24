#!/bin/bash -l
#SBATCH --job-name=tin_expA
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=12:00:00
#SBATCH --output=tin_expA_%j.out
# Add cluster-specific lines (e.g. #SBATCH --account=...) as needed.
#
# Transcript Integrity Number (RSeQC tin.py) for the 7 ExpA replicates, from the paired-end BAMs of
# 00c_trim_star_featurecounts.sh. Run from the repository root:
#   sbatch 08_hidden_degradation/scripts/08a_tin_cluster.sh      (or: bash ... on any machine)
# Needs: gffread (GFF3 -> BED12), samtools, RSeQC (tin.py). Module names are site-specific: set
#   SAMTOOLS_MODULE, GFFREAD_MODULE, PYTHON_MODULE (empty = tool already on PATH) and RSEQC_VENV (a Python
#   virtual environment with RSeQC installed: python3 -m venv <dir>; <dir>/bin/pip install RSeQC).
# EXPERIMENTS selects the experiments (default ExpA; topics 09 and 10 also need ExpB:
#   EXPERIMENTS="ExpA ExpB" sbatch --export=ALL 08_hidden_degradation/scripts/08a_tin_cluster.sh).
# Output: pipeline_data/TIN/<run>_Aligned.sortedByCoord.out.summary.txt and .tin.xls (per transcript), and
#   pipeline_data/TIN/TIN_summary.csv. Skips samples that are already done.
set -euo pipefail
source "$(dirname "$0")/../../config.sh"

BAM_DIR="${BAM_DIR:-$PIPELINE_DIR/STAR_both}"
TIN_DIR="${TIN_DIR:-$PIPELINE_DIR/TIN}"
ANNOT_BED12="${ANNOT_GFF3%.gff3}.bed12"
MAP="data/replicate_map.tsv"
mkdir -p "$TIN_DIR"

load() { [ -n "${1:-}" ] && module add "$1" || true; }
load "${GFFREAD_MODULE:-}"
load "${SAMTOOLS_MODULE:-}"
load "${PYTHON_MODULE:-}"
[ -n "${RSEQC_VENV:-}" ] && source "$RSEQC_VENV/bin/activate"
for tool in gffread samtools tin.py; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found" >&2; exit 1; }
done

# GFF3 -> BED12 (once)
if [ ! -s "$ANNOT_BED12" ]; then
  gffread "$ANNOT_GFF3" --bed -o "$ANNOT_BED12"
  [ -s "$ANNOT_BED12" ] || { echo "ERROR: empty BED12" >&2; exit 1; }
fi

# runs of the selected experiments (ExpA = paper replicates 1-7)
EXPERIMENTS="${EXPERIMENTS:-ExpA}"
RUNS=$(awk -F'\t' -v e=" $EXPERIMENTS " 'NR > 1 && index(e, " " $3 " ") {print $1}' "$MAP")
for run in $RUNS; do
  bam="$BAM_DIR/${run}_Aligned.sortedByCoord.out.bam"
  [ -f "$bam" ] || { echo "missing $bam, skipping" >&2; continue; }
  [ -f "$TIN_DIR/${run}_Aligned.sortedByCoord.out.summary.txt" ] && continue
  [ -f "$bam.bai" ] || samtools index -@ "$THREADS" "$bam"
  ( cd "$TIN_DIR" && tin.py -i "$(realpath "$bam")" -r "$(realpath "$ANNOT_BED12")" )   # tin.py writes to the working directory
done

# one table: run, mean, median and SD of transcript TIN
out="$TIN_DIR/TIN_summary.csv"
echo "run,TIN_mean,TIN_median,TIN_stdev" > "$out"
for run in $RUNS; do
  f="$TIN_DIR/${run}_Aligned.sortedByCoord.out.summary.txt"
  [ -f "$f" ] && tail -n +2 "$f" | awk -v r="$run" -F'\t' 'BEGIN{OFS=","} {print r, $2, $3, $4}' >> "$out"
done
cat "$out"
