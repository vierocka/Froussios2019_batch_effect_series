#!/bin/bash
# Align raw reads with the original paper's STAR settings and count reads on the annotated rRNA genes.
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
ulimit -n 65536 2>/dev/null || true

OUT="$WORK_DIR/rRNA_qc"
STAR_DIR="$OUT/STAR_raw"
STAR_INDEX_PAPER="$REF_DIR/STAR_index_TAIR10_sjdb99"
RRNA_GFF3="$REF_DIR/Arabidopsis_thaliana.TAIR10.63.rRNA_only.gff3"
mkdir -p "$STAR_DIR" "$OUT/counts"

# rRNA loci are typed ncRNA_gene in this annotation
[ -s "$RRNA_GFF3" ] || grep ribos "$ANNOT_GFF3" | grep ncRNA_gene > "$RRNA_GFF3"
[ -s "$RRNA_GFF3" ] || { echo "no rRNA genes found in annotation" >&2; exit 1; }

if [ ! -f "$STAR_INDEX_PAPER/SAindex" ]; then
  mkdir -p "$STAR_INDEX_PAPER"
  STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$STAR_INDEX_PAPER" \
       --genomeFastaFiles "$GENOME_FA" --sjdbGTFfile "$ANNOT_GFF3" --sjdbOverhang 99 --genomeSAindexNbases 12
fi

for run in $(tail -n +2 "$RUN_LIST" | cut -f1); do
  R1="$FASTQ_DIR/${run}_1.fastq.gz"; R2="$FASTQ_DIR/${run}_2.fastq.gz"
  [ -f "$R1" ] && [ -f "$R2" ] || { echo "missing fastq for $run, skipping" >&2; continue; }
  BAM="$STAR_DIR/${run}_Aligned.sortedByCoord.out.bam"

  [ -f "$BAM" ] || STAR --runThreadN "$THREADS" --genomeDir "$STAR_INDEX_PAPER" \
    --outSAMstrandField intronMotif --outSJfilterIntronMaxVsReadN 5000 10000 15000 20000 \
    --outFilterType BySJout --outFilterMultimapNmax 2 --outFilterMismatchNmax 5 \
    --outFileNamePrefix "$STAR_DIR/${run}_" --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 100000000000 --readFilesCommand gunzip -c --readFilesIn "$R1" "$R2"

  [ -f "$OUT/counts/${run}_rRNA.count" ] || featureCounts -T "$THREADS" -t ncRNA_gene -g gene_id -s 2 -p -P -B \
    -a "$RRNA_GFF3" -o "$OUT/counts/${run}_rRNA.count" "$BAM"
done
