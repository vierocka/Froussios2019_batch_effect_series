#!/bin/bash
# Download the 17 runs of ERP021226 with sra-tools. Skips runs already done.
# Usage: bash 00_download/download_sra.sh [RUN ...]   (default: all runs in data/run_list.tsv)
set -euo pipefail
source "$(dirname "$0")/../../config.sh"

mkdir -p "$FASTQ_DIR" "$WORK_DIR/sra_tmp"

if [ "$#" -gt 0 ]; then RUNS=("$@"); else RUNS=($(tail -n +2 "$RUN_LIST" | cut -f1)); fi

for run in "${RUNS[@]}"; do
  if [ -s "$FASTQ_DIR/${run}_1.fastq.gz" ] && [ -s "$FASTQ_DIR/${run}_2.fastq.gz" ]; then
    echo "skip $run"; continue
  fi
  prefetch --max-size 100G -O "$WORK_DIR/sra_tmp" "$run"
  fasterq-dump --split-files --threads "$THREADS" -O "$WORK_DIR/sra_tmp" -t "$WORK_DIR/sra_tmp" "$WORK_DIR/sra_tmp/$run"

  # mates must match
  n1=$(wc -l < "$WORK_DIR/sra_tmp/${run}_1.fastq"); n2=$(wc -l < "$WORK_DIR/sra_tmp/${run}_2.fastq")
  [ "$n1" -eq "$n2" ] || { echo "ERROR: $run mate counts differ ($n1 vs $n2)" >&2; exit 1; }

  gzip -c "$WORK_DIR/sra_tmp/${run}_1.fastq" > "$FASTQ_DIR/${run}_1.fastq.gz"
  gzip -c "$WORK_DIR/sra_tmp/${run}_2.fastq" > "$FASTQ_DIR/${run}_2.fastq.gz"
  rm -rf "$WORK_DIR/sra_tmp/${run}" "$WORK_DIR/sra_tmp/${run}"_[12].fastq
done
