#!/bin/bash

# Parse command-line arguments
while getopts d:r:f: flag
do
    case "${flag}" in
        d) directory=${OPTARG};;
        r) reference_directory=${OPTARG};;
        f) merged_bed_file=${OPTARG};;  # Allow file input
    esac
done

# Check required parameters
if [ -z "$directory" ]; then
    echo "Error: Path to sample directories not provided [-d]"
    exit 1
fi

if [ -z "$reference_directory" ]; then
    echo "Error: Path to reference directory not provided [-r]"
    exit 1
fi

directory=$(echo "${directory}" | sed 's![^/]$!&/!')
reference_directory=$(echo "${reference_directory}" | sed 's![^/]$!&/!')


merged_bed_file="${directory}peaks.family.bed"
output_file="${directory}shared_targets.csv"

if [ ! -f "$merged_bed_file" ]; then
    echo "Error: Merged BED file '$merged_bed_file' not found!"
    exit 1
fi

location=$(conda info | awk '/base environment/' | awk '{print $4}')
source ${location}/etc/profile.d/conda.sh
conda activate SCRAP

# Find overlapping regions within the merged BED file (self-intersection)
echo "Finding shared miRNA & isomiR targets..."
bedtools intersect -a "$merged_bed_file" -b "$merged_bed_file" -wo | awk '$4 != $10' > "${directory}shared_targets_raw.bed"

# Check if any shared targets were found
if [ ! -s "${directory}shared_targets_raw.bed" ]; then
    echo "No shared targets found!"
    exit 0
fi

# Format output to CSV (extracts relevant columns)
echo "Processing shared targets..."
awk 'BEGIN{OFS=","} {print $1, $2, $3, $4, $6, $10, $11, $12}' "${directory}shared_targets_raw.bed" > "$output_file"

# Clean up temporary files
rm "${directory}shared_targets_raw.bed"

echo "Shared miRNA & isomiR targets saved in: $output_file"