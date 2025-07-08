#!/bin/bash

# Usage: ./create_bam.sh <case_number>

CASE_NUM=$1
GENOME_INDEX="./genome/unis/uni"
RAW_DATA_DIR="./raw_data"
BAM_DIR="./bam"

echo "Starting BAM creation for case ${CASE_NUM}..."

# Mapping father, adding @RG tag
echo "Aligning father reads..."
bowtie2 -x $GENOME_INDEX -U $RAW_DATA_DIR/case${CASE_NUM}_father.fq.gz \
--rg-id father${CASE_NUM} --rg "SM:father${CASE_NUM}" --rg "PL:ILLUMINA" | \
samtools view -bSu - | samtools sort -@ 4 -o $BAM_DIR/case${CASE_NUM}_father.bam
samtools index -@ 4 $BAM_DIR/case${CASE_NUM}_father.bam
echo "Father mapped, sorted, and indexed."

# Mapping mother, adding @RG tag
echo "Aligning mother reads..."
bowtie2 -x $GENOME_INDEX -U $RAW_DATA_DIR/case${CASE_NUM}_mother.fq.gz \
--rg-id mother${CASE_NUM} --rg "SM:mother${CASE_NUM}" --rg "PL:ILLUMINA" | \
samtools view -bSu - | samtools sort -@ 4 -o $BAM_DIR/case${CASE_NUM}_mother.bam
samtools index -@ 4 $BAM_DIR/case${CASE_NUM}_mother.bam
echo "Mother mapped, sorted, and indexed."

# Mapping child, adding @RG tag
echo "Aligning child reads..."
bowtie2 -x $GENOME_INDEX -U $RAW_DATA_DIR/case${CASE_NUM}_child.fq.gz \
--rg-id child${CASE_NUM} --rg "SM:child${CASE_NUM}" --rg "PL:ILLUMINA" | \
samtools view -bSu - | samtools sort -@ 4 -o $BAM_DIR/case${CASE_NUM}_child.bam
samtools index -@ 4 $BAM_DIR/case${CASE_NUM}_child.bam
echo "Child mapped, sorted, and indexed."

echo "BAM files for case ${CASE_NUM} created and indexed."
