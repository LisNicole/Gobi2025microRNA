#!/bin/sh

# Parse arguments
while getopts d:a:p: flag
do
    case "${flag}" in
        d) directory=${OPTARG};;  # Sample directory
        a) adapter_file=${OPTARG};;  # Adapter file
        p) paired_end=${OPTARG};;  # Paired-end sequencing (yes/no)
    esac
done

# Check if all required arguments are provided
if [ -z "$directory" ]; then
    echo "Error: Path to sample directories not provided [-d]"
    exit 1
fi
if [ -z "$adapter_file" ]; then
    echo "Error: Path to adapter file not provided [-a]"
    exit 1
fi
if [ -z "$paired_end" ]; then
    echo "Error: Are samples paired-end? [-p]"
    exit 1
fi

# Ensure paths end with /
directory=$(echo "${directory}" | sed 's![^/]$!&/!')

# Get sample name from adapter file
sample=$(awk '!/^#/ {print $1; exit}' ${adapter_file})

echo "Running test on sample: ${sample}"

### Download GtRNAdb.fasta ###
echo "Downloading GtRNAdb.fasta file"
wget https://gtrnadb.ucsc.edu/download/gtRNAdb.fasta -O ${directory}GtRNAdb.fasta

echo "GtRNAdb.fasta downloaded successfully!"

echo "Test completed successfully!"


for sample in $samples
do

  # Extract adapter and barcode sequences from adapter file
  five_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $2}')
  three_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $3}')
  five_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $4}')
  three_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $5}')
