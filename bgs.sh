CASE_NUM=$1

bedtools genomecov -ibam ./bam/case${CASE_NUM}_father.bam -bg \
-trackline -trackopts 'name="father"' -max 100 > \
./bgs/case${CASE_NUM}_fatherCov.bg


bedtools genomecov -ibam ./bam/case${CASE_NUM}_mother.bam -bg \
-trackline -trackopts 'name="mother"' -max 100 > \
./bgs/case${CASE_NUM}_motherCov.bg


bedtools genomecov -ibam ./bam/case${CASE_NUM}_child.bam -bg \
-trackline -trackopts 'name="child"' -max 100 > \
./bgs/case${CASE_NUM}_childCov.bg
