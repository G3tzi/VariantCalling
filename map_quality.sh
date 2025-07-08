#!/bin/bash

# Usage: ./map_quality.sh <case_number>

CASE_NUM=$1
BAM_DIR="./bam"
GENOME_DIR="./genome"
MAP_QUALITY_DIR="./map_quality"

echo "Starting map quality assessment for case ${CASE_NUM}..."

# Perform Qualimap for father
echo "Running Qualimap for father..."
qualimap bamqc -bam $BAM_DIR/case${CASE_NUM}_father.bam -gff $GENOME_DIR/*.bed -outdir $MAP_QUALITY_DIR/case${CASE_NUM}_father_qualimap

# Perform Qualimap for mother
echo "Running Qualimap for mother..."
qualimap bamqc -bam $BAM_DIR/case${CASE_NUM}_mother.bam -gff $GENOME_DIR/*.bed -outdir $MAP_QUALITY_DIR/case${CASE_NUM}_mother_qualimap

# Perform Qualimap for child
echo "Running Qualimap for child..."
qualimap bamqc -bam $BAM_DIR/case${CASE_NUM}_child.bam -gff $GENOME_DIR/*.bed -outdir $MAP_QUALITY_DIR/case${CASE_NUM}_child_qualimap

echo "Map quality metrics calculated for all samples."
