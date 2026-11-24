#!/bin/bash
# Recount the existing STAR_both BAMs with featureCounts -s 0, -s 1 and -s 2 (and -s 2 counting
# fragments instead of reads). No realignment. Flags otherwise identical to 00c_trim_star_featurecounts.sh,
# so -s is the only thing that changes between fC_s0, fC_s1 and fC_s2.
#   fC_s0         unstranded; identical to fC/
#   fC_s1         stranded, read 1 sense
#   fC_s2         stranded, read 1 antisense (expected for TruSeq Stranded Total RNA)
#   fC_s2_pairs   as fC_s2, plus --countReadPairs: one count per fragment, not per mate
#                 (featureCounts >= 2.0.2 counts each mate separately with -p alone)
# Skips outputs that already exist. Usage: bash 00_preprocessing/scripts/00e_count_strandedness.sh
set -euo pipefail
source "$(dirname "$0")/../../config.sh"

BAM_DIR="${BAM_DIR:-$PIPELINE_DIR/STAR_both}"
OUT_BASE="${OUT_BASE:-$PIPELINE_DIR/STAR_both}"

declare -A OPTS=( [fC_s0]="-s 0" [fC_s1]="-s 1" [fC_s2]="-s 2" [fC_s2_pairs]="-s 2 --countReadPairs" )

for run in $(tail -n +2 "$RUN_LIST" | cut -f1); do
  BAM="$BAM_DIR/${run}_Aligned.sortedByCoord.out.bam"
  [ -f "$BAM" ] || { echo "missing $BAM, skipping" >&2; continue; }
  for v in fC_s0 fC_s1 fC_s2 fC_s2_pairs; do
    mkdir -p "$OUT_BASE/$v"
    out="$OUT_BASE/$v/${run}.count"
    [ -f "$out" ] && [ -f "$out.summary" ] && continue
    # shellcheck disable=SC2086
    featureCounts -T "$THREADS" -p ${OPTS[$v]} -a "$ANNOT_GFF3" -t gene -g gene_id -o "$out" "$BAM"
  done
done

# one-line overview: assigned fraction per run and variant
for v in fC_s0 fC_s1 fC_s2 fC_s2_pairs; do
  for s in "$OUT_BASE/$v"/*.count.summary; do
    awk -v v="$v" -v f="$(basename "$s" .count.summary)" 'NR>1{t+=$2} $1=="Assigned"{a=$2}
      END{if (t > 0) printf "%s\t%s\tassigned %d of %d (%.1f%%)\n", v, f, a, t, 100*a/t}' "$s"
  done
done
