#!/bin/bash

CASE_NUM=$1
VCF_DIR="./vcfs"

samples=$(bcftools query -l ${VCF_DIR}/case${CASE_NUM}_variants.vcf)

sample_array=($samples)

mother_index=-1
father_index=-1
child_index=-1

for i in "${!sample_array[@]}"; do
	if [[ "${sample_array[$i]}" == "mother${CASE_NUM}" ]]; then
		mother_index=$i
	elif [[ "${sample_array[$i]}" == "father${CASE_NUM}" ]]; then
		father_index=$i
	elif [[ "${sample_array[$i]}" == "child${CASE_NUM}" ]]; then
		child_index=$i
	fi
done

if [[ $mother_index -eq -1 || $father_index -eq -1 || $child_index -eq -1 ]]; then
	echo "Error: Could not find the required samples (mother${CASE_NUM}, father${CASE_NUM}, child${CASE_NUM}"
	exit 1
fi

FILTER="GT[$mother_index]='0/0' && GT[$father_index]='0/0' && GT[$child_index]='0/1'"

bcftools view -i "$FILTER" ${VCF_DIR}/case${CASE_NUM}_targeted.vcf -o ${VCF_DIR}/case${CASE_NUM}_dominant.vcf

FILTER="GT[$mother_index]='0/1' && GT[$father_index]='0/1' && GT[$child_index]='1/1'"

bcftools view -i "$FILTER" ${VCF_DIR}/case${CASE_NUM}_targeted.vcf -o ${VCF_DIR}/case${CASE_NUM}_recessive.vcf

echo "Filtered vcf for case: ${CASE_NUM}"
