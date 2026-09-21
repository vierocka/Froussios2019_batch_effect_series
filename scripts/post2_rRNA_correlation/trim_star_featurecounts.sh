#!/bin/bash
# Trimmomatic -> STAR (both-mates-survived pairs, and forward-only singletons) -> featureCounts.
# Skips every step whose output already exists.
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
ulimit -n 65536 2>/dev/null || true

ADAPTERS="${ADAPTERS:?set ADAPTERS to the path of TruSeq3-PE.fa}"
TRIM_DIR="$WORK_DIR/trimmomatic"
PAIRED_DIR="$WORK_DIR/STAR_both"
FONLY_DIR="$WORK_DIR/STAR_Fonly"
mkdir -p "$TRIM_DIR" "$PAIRED_DIR/fC" "$FONLY_DIR/fC"

STAR_OPTS=(--runThreadN "$THREADS" --genomeDir "$STAR_INDEX"
  --outFilterType BySJout --outFilterMultimapNmax 10
  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
  --outFilterMismatchNmax 12 --outFilterMismatchNoverReadLmax 0.04
  --alignIntronMin 20 --outSAMtype BAM SortedByCoordinate
  --limitBAMsortRAM 12000000000 --outBAMsortingBinsN 20 --readFilesCommand gunzip -c)

for run in $(tail -n +2 "$RUN_LIST" | cut -f1); do
  R1="$FASTQ_DIR/${run}_1.fastq.gz"; R2="$FASTQ_DIR/${run}_2.fastq.gz"
  [ -f "$R1" ] && [ -f "$R2" ] || { echo "missing fastq for $run, skipping" >&2; continue; }

  if [ ! -f "$TRIM_DIR/${run}_1P.fq.gz" ]; then
    trimmomatic PE -threads "$THREADS" -summary "$TRIM_DIR/${run}_summary.txt" -quiet -validatePairs \
      "$R1" "$R2" \
      "$TRIM_DIR/${run}_1P.fq.gz" "$TRIM_DIR/${run}_1U.fq.gz" \
      "$TRIM_DIR/${run}_2P.fq.gz" "$TRIM_DIR/${run}_2U.fq.gz" \
      ILLUMINACLIP:"$ADAPTERS":2:28:8 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
  fi

  BAM_P="$PAIRED_DIR/${run}_Aligned.sortedByCoord.out.bam"
  BAM_F="$FONLY_DIR/${run}_Aligned.sortedByCoord.out.bam"
  [ -f "$BAM_P" ] || STAR "${STAR_OPTS[@]}" --outFileNamePrefix "$PAIRED_DIR/${run}_" \
    --readFilesIn "$TRIM_DIR/${run}_1P.fq.gz" "$TRIM_DIR/${run}_2P.fq.gz"
  [ -f "$BAM_F" ] || STAR "${STAR_OPTS[@]}" --outFileNamePrefix "$FONLY_DIR/${run}_" \
    --readFilesIn "$TRIM_DIR/${run}_1U.fq.gz"

  [ -f "$PAIRED_DIR/fC/${run}.count" ] || featureCounts -T "$THREADS" -p -a "$ANNOT_GFF3" -t gene -g gene_id \
    -o "$PAIRED_DIR/fC/${run}.count" "$BAM_P"
  [ -f "$FONLY_DIR/fC/${run}.count" ]  || featureCounts -T "$THREADS" -a "$ANNOT_GFF3" -t gene -g gene_id \
    -o "$FONLY_DIR/fC/${run}.count" "$BAM_F"
done
