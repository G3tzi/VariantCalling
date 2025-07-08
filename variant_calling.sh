#!/bin/bash

# Usage: ./variant_calling.sh <case_number>

CASE_NUM=$1
BAM_DIR="./bam"
VCF_DIR="./vcfs"
GENOME_FASTA="./genome/universe.fasta"
EXONS_BED="./genome/exons16Padded_sorted.bed"

echo "Starting variant calling for case ${CASE_NUM}..."

# Run FreeBayes for variant calling
echo "Calling variants using FreeBayes..."
freebayes -f $GENOME_FASTA \
-m 20 -C 5 -Q 10 --min-coverage 10 \
$BAM_DIR/case${CASE_NUM}_father.bam $BAM_DIR/case${CASE_NUM}_mother.bam $BAM_DIR/case${CASE_NUM}_child.bam \
-t $EXONS_BED \
> $VCF_DIR/case${CASE_NUM}_variants.vcf

echo "Variant calling complete. VCF file created: ${VCF_DIR}/${CASE_NUM}_variants.vcf"
