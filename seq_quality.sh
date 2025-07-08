CASE_NUM=$1
RAW_DATA_DIR="./raw_data/"
SEQ_QUALITY_DIR="./seq_quality"

mkdir -p $SEQ_QUALITY_DIR

echo "Starting sequencing quality analysis for case ${CASE_NUM}..."

fastqc -o $SEQ_QUALITY_DIR $RAW_DATA_DIR/case${CASE_NUM}_father.fq.gz
fastqc -o $SEQ_QUALITY_DIR $RAW_DATA_DIR/case${CASE_NUM}_mother.fq.gz
fastqc -o $SEQ_QUALITY_DIR $RAW_DATA_DIR/case${CASE_NUM}_child.fq.gz

echo "Sequencing quality analysis completed. Result stored in $SEQ_QUALITY_DIR."
