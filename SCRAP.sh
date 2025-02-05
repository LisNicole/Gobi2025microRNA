#!/bin/bash

while getopts d:a:p:f:r:m:g: flag
do
    case "${flag}" in
        d) directory=${OPTARG};;  # Sample directory
        a) adapter_file=${OPTARG};;  # Adapter file
        p) paired_end=${OPTARG};;  # Paired-end sequencing (yes/no)
        f) pre_filtered=${OPTARG};;  # Filter pre-miRNAs (yes/no)
        m) miRBase_species_abbreviation=${OPTARG};;  # miRBase species
        g) genome_species_abbreviation=${OPTARG};;  # Genome species
    esac
done

# Check if all required arguments are provided
if [ -z "$directory" ]; then echo "Error: Path to sample directories not provided [-d]"; exit 1; fi
if [ -z "$adapter_file" ]; then echo "Error: Path to adapter file not provided [-a]"; exit 1; fi
if [ -z "$paired_end" ]; then echo "Error: Are samples paired-end? [-p]"; exit 1; fi
if [ -z "$pre_filtered" ]; then echo "Error: Filter out pre-miRNAs and tRNAs? [-f]"; exit 1; fi
if [ -z "$miRBase_species_abbreviation" ]; then echo "Error: miRBase species abbreviation not provided [-m]"; exit 1; fi
if [ -z "$genome_species_abbreviation" ]; then echo "Error: Species genome abbreviation not provided [-g]"; exit 1; fi

# Ensure paths end with /
directory=$(echo "${directory}" | sed 's![^/]$!&/!')
reference_directory=$(echo "${reference_directory}" | sed 's![^/]$!&/!')

# Activate the virtual environment
source .venv/bin/activate

# Get sample name from adapter file (assuming first column contains sample names)
sample=$(awk '!/^#/' ${adapter_file} | awk '{print $1}' | head -n 1)

echo "Running test on sample: ${sample}"

### Quality Control (FastQC) ###
mkdir -p ${directory}FastQC_Reports

for file in $(ls ${directory}${sample}.fastq | awk '/fastq/')

do
    echo "Running FastQC on ${file}"
    fastqc -o ${directory}FastQC_Reports ${directory}${sample}/${file}
done

### Quality Control (MultiQC) ###
mkdir -p ${directory}MultiQC_Report

echo "Running MultiQC"
multiqc ${directory}FastQC_Reports -o ${directory}MultiQC_Report

echo "Test completed successfully!"
