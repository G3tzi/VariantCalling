# Genomics Project: Molecular Diagnosis of Rare Genetic Disorders

## Project Overview

This project aims to diagnose rare Mendelian disorders in simulated TRIO data (father, mother, child) using exome sequencing data from chromosome 16 (hg19). The analysis focuses on identifying disease-causing variants in children who may be affected by autosomal dominant (AD) or recessive (AR) disorders, while parents remain healthy.

### Key Objectives
- Analyze 10 simulated TRIO cases using exome sequencing data
- Identify disease-causing variants following Mendelian inheritance patterns
- Distinguish between autosomal dominant and recessive inheritance modes
- Provide molecular diagnosis for rare genetic disorders

## Dataset

- **Reference Genome**: GRCh37/hg19
- **Target Region**: Chromosome 16 exons only
- **Data Format**: FASTQ files for each trio member (father, mother, child)
- **Inheritance Models**: 
  - Autosomal Dominant (AD): De novo mutations (0/1 in child, 0/0 in parents)
  - Autosomal Recessive (AR): Compound heterozygous (1/1 in child, 0/1 in parents)

## Pipeline Architecture

The analysis pipeline consists of several automated bash scripts that process the data sequentially:

### 1. Data Preparation
- **`link_cases.sh`**: Creates symbolic links to raw FASTQ files
- **Usage**: `./link_cases.sh <case_number>`

### 2. Sequence Alignment
- **`create_bam.sh`**: Aligns reads to reference genome and creates BAM files
- Uses Bowtie2 for alignment with proper read group tags
- Sorts and indexes BAM files using samtools
- **Usage**: `./create_bam.sh <case_number>`

### 3. Quality Control
- **`seq_quality.sh`**: Performs sequencing quality analysis using FastQC
- **`map_quality.sh`**: Assesses mapping quality using Qualimap
- **Usage**: `./seq_quality.sh <case_number>` or `./map_quality.sh <case_number>`

### 4. Variant Calling
- **`variant_calling.sh`**: Calls variants using FreeBayes
- Applies quality filters and targets exonic regions
- **Usage**: `./variant_calling.sh <case_number>`

### 5. Variant Filtering and Annotation
- **`annotations.sh`**: Filters variants by quality and intersects with target regions
- **`create_recessive_dominant_vcfs.sh`**: Separates variants by inheritance pattern
- **Usage**: `./annotations.sh <case_number>` and `./create_recessive_dominant_vcfs.sh <case_number>`

### 6. Coverage Analysis
- **`bgs.sh`**: Generates coverage tracks for UCSC Genome Browser visualization
- **Usage**: `./bgs.sh <case_number>`

### 7. Master Pipeline
- **`run_pipeline.sh`**: Executes the complete analysis pipeline
- **Usage**: `./run_pipeline.sh <case_number> [n]`
  - Add `n` as second parameter to run without interactive prompts

## Key Analysis Steps

### Variant Prioritization Strategy

1. **Quality Filtering**: 
   - Variant quality (QUAL) > 30
   - Sample depth (FMT/DP) ≥ 15
   - Overall depth (INFO/DP) ≥ 20

2. **Target Region Intersection**: 
   - Variants retained only if within exonic regions
   - Uses BEDTools intersection with `exons16Padded_sorted.bed`

3. **Inheritance Pattern Filtering**:
   - **Autosomal Dominant**: `GT[mother]='0/0' && GT[father]='0/0' && GT[child]='0/1'`
   - **Autosomal Recessive**: `GT[mother]='0/1' && GT[father]='0/1' && GT[child]='1/1'`

### Annotation and Visualization

- **VEP Annotation**: Functional annotation using Variant Effect Predictor
- **UCSC Genome Browser**: Visualization of coverage tracks and variant locations
- **gnomAD Filtering**: Allele frequency ≤ 1×10⁻⁴ for rare variants

## Results Summary

The pipeline successfully identified disease-causing variants in 10 cases:

| Cases | Gene | Inheritance | Disease |
|-------|------|-------------|---------|
| 617, 628, 596 | FANCA | AR | Fanconi Anemia |
| 657 | CYLD | AD | Familial cylindromatosis |
| 717 | ANKRD11 | AD | KBG syndrome |
| 737, 639, 642, 645, 701 | CREBBP | AD | Rubinstein-Taybi syndrome |

### Variant Types Identified:
- Frameshift mutations
- Stop-gained variants
- Splice acceptor variants

## How to Run the Pipeline

### Prerequisites
- Access to Unix/Linux environment with required tools installed:
  - Bowtie2, samtools, bcftools
  - FreeBayes, FastQC, Qualimap
  - BEDTools, MultiQC

### Quick Start

1. **Set up directory structure**:
```bash
mkdir -p raw_data bam vcfs seq_quality map_quality multiqc_reports bgs
```

2. **Link raw data**:
```bash
./link_cases.sh <case_number>
```

3. **Run complete pipeline**:
```bash
./run_pipeline.sh <case_number>
```

4. **Or run individual steps**:
```bash
./create_bam.sh <case_number>
./map_quality.sh <case_number>
./variant_calling.sh <case_number>
./seq_quality.sh <case_number>
./annotations.sh <case_number>
./create_recessive_dominant_vcfs.sh <case_number>
./bgs.sh <case_number>
```

### Output Files

- **BAM files**: `./bam/case<N>_*.bam`
- **VCF files**: `./vcfs/case<N>_*.vcf`
- **Quality reports**: `./multiqc_reports/case<N>/`
- **Coverage tracks**: `./bgs/case<N>_*Cov.bg`

## Technical Implementation

### Key Tools and Parameters

- **Bowtie2**: `--rg-id` and `--rg` for read group assignment
- **FreeBayes**: `-m 20 -C 5 -Q 10 --min-coverage 10` for variant calling
- **BCFTools**: Quality filtering and genotype-based inheritance filtering
- **BEDTools**: Genomic interval operations for target region intersection

### Quality Control Metrics

- High alignment rates (>99%)
- Sufficient mean coverage across samples
- GC content and duplication levels within expected ranges
- Phred scores >30 for reliable base calling

## Project Structure

```
genomics_project/
├── raw_data/          # Symbolic links to FASTQ files
├── bam/               # Aligned BAM files
├── vcfs/              # Variant call files
├── seq_quality/       # FastQC reports
├── map_quality/       # Qualimap reports
├── multiqc_reports/   # Aggregated QC reports
├── bgs/               # Coverage tracks
└── *.sh               # Pipeline scripts
```

## Authors

- Gabriele Ghezzi
- Lorenzo Paratici

**Course**: Genomics 2025, University of Milan  
**Instructor**: Prof. Matteo Chiara  
**Date**: May 4, 2025
