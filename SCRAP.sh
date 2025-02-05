#!/bin/bash

while getopts d:a:p: flag
do
    case "${flag}" in
        d) directory=${OPTARG};;  # Sample directory
        a) adapter_file=${OPTARG};;  # Adapter file
        p) paired_end=${OPTARG};;  # Paired-end sequencing (yes/no)
    esac
done

# Check if all required arguments are provided
if [ -z "$directory" ]; then echo "Error: Path to sample directories not provided [-d]"; exit 1; fi
if [ -z "$adapter_file" ]; then echo "Error: Path to adapter file not provided [-a]"; exit 1; fi
if [ -z "$paired_end" ]; then echo "Error: Are samples paired-end? [-p]"; exit 1; fi

# Ensure paths end with /
directory=$(echo "${directory}" | sed 's![^/]$!&/!')

# Activate the virtual environment
source .venv/bin/activate

# Get sample name from adapter file
sample=$(awk '!/^#/ {print $1; exit}' ${adapter_file})

echo "Running test on sample: ${sample}"

### Quality Control (FastQC) ###
mkdir -p ${directory}FastQC_Reports

for file in ${directory}${sample}*.fastq
do
    echo "Running FastQC on ${file}"
    fastqc -o ${directory}FastQC_Reports "${file}"
done

### Quality Control (MultiQC) ###
mkdir -p ${directory}MultiQC_Report

echo "Running MultiQC"
multiqc ${directory}FastQC_Reports -o ${directory}MultiQC_Report

echo "Test completed successfully!"
