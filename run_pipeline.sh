#!/bin/bash

# Usage: ./run_pipeline.sh <case_number>

CASE_NUM=$1
CHECK=$2

if [ -z "$CASE_NUM" ]; then
    echo "Error: Please provide a case number (e.g., ./run_pipeline.sh <case_number>)"
    exit 1
fi

echo "Starting pipeline for case ${CASE_NUM}..."

# Run the create_bam script
echo "Step 1: Creating BAM files..."
./create_bam.sh $CASE_NUM
echo "BAM files created for case ${CASE_NUM}."

if [["$CHECK" != "n"]]; then
	read -p "Do you want to proceed to map quality assessment? (y/n):  " proceed
	if [["$proceed" != "y"]]; then
		echo "Exiting pipeline."
		exit 0
	fi
fi

# Run the map_quality script
echo "Step 2: Running map quality assessment..."
./map_quality.sh $CASE_NUM
echo "Map quality assessment complete for case ${CASE_NUM}."

if [["$CHECK" != "n"]]; then
	read -p "Do you want to proceed to variant calling? (y/n): " proceed
	if [["$proceed" != "y"]]; then
		echo "Exiting pipeline."
		exit 0
	fi
fi

# Run the variant_calling script
echo "Step 3: Performing variant calling..."
./variant_calling.sh $CASE_NUM
echo "Variant calling complete for case ${CASE_NUM}."

# Run the seq_quality script
echo "Step 4: Sequencing Quality Analysis"
./seq_quality.sh $CASE_NUM

echo "Step 5: Running MultiQC..."
multiqc ./seq_quality/case${CASE_NUM}* ./map_quality/case${CASE_NUM}* -o ./multiqc_reports/case${CASE_NUM}
echo "MultiQC report created for case ${CASE_NUM}."

./annotations.sh ${CASE_NUM}
./create_recessive_dominant.sh ${CASE_NUM}

echo "Pipeline completed for case ${CASE_NUM}. All steps are done."
