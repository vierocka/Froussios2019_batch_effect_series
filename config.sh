# Shared settings; override any of these via environment variables.
export WORK_DIR="${WORK_DIR:-$PWD/work}"
export THREADS="${THREADS:-8}"
export FASTQ_DIR="${FASTQ_DIR:-$WORK_DIR/fastq}"
export REF_DIR="${REF_DIR:-$WORK_DIR/ref}"
export RUN_LIST="${RUN_LIST:-$PWD/data/run_list.tsv}"
export GENOME_FA="$REF_DIR/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
export ANNOT_GFF3="$REF_DIR/Arabidopsis_thaliana.TAIR10.63.gff3"
export STAR_INDEX="$REF_DIR/STAR_index_TAIR10"
