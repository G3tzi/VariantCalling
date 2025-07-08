CASE_NUM=$1
ln -s /home/BCG2025_genomics_exam/case${CASE_NUM}_child.fq.gz ./raw_data/case${CASE_NUM}_child.fq.gz
ln -s /home/BCG2025_genomics_exam/case${CASE_NUM}_father.fq.gz ./raw_data/case${CASE_NUM}_father.fq.gz
ln -s /home/BCG2025_genomics_exam/case${CASE_NUM}_mother.fq.gz ./raw_data/case${CASE_NUM}_mother.fq.gz
