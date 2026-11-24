# Shared settings; override any of these via environment variables.
export PIPELINE_DIR="${PIPELINE_DIR:-$PWD/pipeline_data}"
export THREADS="${THREADS:-8}"
export FASTQ_DIR="${FASTQ_DIR:-$PIPELINE_DIR/fastq}"
export REF_DIR="${REF_DIR:-$PIPELINE_DIR/ref}"
export RUN_LIST="${RUN_LIST:-$PWD/data/run_list.tsv}"
export GENOME_FA="$REF_DIR/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
export ANNOT_GFF3="$REF_DIR/Arabidopsis_thaliana.TAIR10.63.gff3"
export STAR_INDEX="$REF_DIR/STAR_index_TAIR10"
