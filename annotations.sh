
CASE_NUM=$1

bcftools filter -i 'QUAL > 30 && FMT/DP>=15 && INFO/DP>=20' ./vcfs/case${CASE_NUM}_variants.vcf | \
bedtools intersect -a stdin -b ./genome/exons16Padded_sorted.bed -header -u > ./vcfs/case${CASE_NUM}_targeted.vcf

