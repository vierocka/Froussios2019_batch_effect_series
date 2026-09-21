#!/bin/bash
# Download Ensembl Plants release 63 TAIR10 genome + annotation and build the STAR index.
set -euo pipefail
source "$(dirname "$0")/../../config.sh"
mkdir -p "$REF_DIR"

BASE=https://ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-63
[ -f "$GENOME_FA" ]  || { curl -sS -o "$GENOME_FA.gz"  "$BASE/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"; gunzip "$GENOME_FA.gz"; }
[ -f "$ANNOT_GFF3" ] || { curl -sS -o "$ANNOT_GFF3.gz" "$BASE/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3.gz"; gunzip "$ANNOT_GFF3.gz"; }

# gene_id must exist on gene records (featureCounts -g gene_id)
awk -F'\t' '$3=="gene"{print; exit}' "$ANNOT_GFF3" | grep -q gene_id || { echo "no gene_id attribute" >&2; exit 1; }

if [ ! -f "$STAR_INDEX/SAindex" ]; then
  mkdir -p "$STAR_INDEX"
  STAR --runThreadN "$THREADS" --runMode genomeGenerate \
       --genomeDir "$STAR_INDEX" --genomeFastaFiles "$GENOME_FA" \
       --sjdbGTFfile "$ANNOT_GFF3" --sjdbGTFtagExonParentTranscript Parent \
       --sjdbGTFfeatureExon exon --sjdbOverhang 100 --genomeSAindexNbases 12
fi
