#!/bin/bash

# Input arguments
while getopts d:r:t:o: flag
do
    case "${flag}" in
        d) directory=${OPTARG};;  # Working directory
        r) reference_utr=${OPTARG};;  # 3'UTR reference BED file
        t) target_file=${OPTARG};;  # miRNA target BED file
        o) output_file=${OPTARG};;  # Output file (Optional)
    esac
done


if [[ -z "$directory" || -z "$reference_utr" || -z "$target_file" ]]; then
    echo "Usage: $0 -d <directory> -r <3UTR_reference> -t <target_file>"
    exit 1
fi


if [[ -z "$output_file" ]]; then
    output_file="$directory/shared_3utr_targets.csv"
fi

echo "Filtering miRNA targets within 3'UTRs..."
bedtools intersect -a "$target_file" -b "$reference_utr" -wa > "$directory/filtered_3UTR_targets.bed"

echo "Grouping miRNAs by family..."
awk '{print $4}' "$directory/filtered_3UTR_targets.bed" | cut -d'|' -f1 | sort | uniq > "$directory/unique_miRNA_families.txt"

echo "Computing Jaccard index for shared targets..."
echo "miRNA,Jaccard_Index" > "$output_file"

while read mir_family; do
    # Get targets for the miRNA family
    mir_targets=$(awk -v mir="$mir_family" '$4 ~ mir {print $1,$2,$3}' "$directory/filtered_3UTR_targets.bed")

    # Find isomiRs (other variants of the same family) in the BED file
    iso_targets=$(awk -v mir="$mir_family" '$4 ~ mir {print $1,$2,$3}' "$directory/filtered_3UTR_targets.bed")

    # Calculate Jaccard index
    shared_targets=$(echo -e "$mir_targets\n$iso_targets" | sort | uniq -d | wc -l)
    total_mir_targets=$(echo "$mir_targets" | wc -l)

    if [[ $total_mir_targets -gt 0 ]]; then
        jaccard_index=$(echo "scale=2; ($shared_targets / $total_mir_targets) * 100" | bc)
        echo "$mir_family,$jaccard_index" >> "$output_file"
    fi
done < "$directory/unique_miRNA_families.txt"

echo "Generating histogram..."
python3 - <<END
import pandas as pd
import matplotlib.pyplot as plt

# Load Jaccard index data
df = pd.read_csv("$output_file")
jaccard_values = df["Jaccard_Index"]

# Plot histogram
plt.hist(jaccard_values, bins=10, edgecolor='black', alpha=0.7)
plt.xlabel("Jaccard Index (%)")
plt.ylabel("Count")
plt.title("Histogram of Shared miRNA-IsomiR Targets (Jaccard Index)")
plt.savefig("$directory/shared_targets_histogram.png")
plt.show()
END

echo "Analysis complete. Results saved in $output_file and shared_targets_histogram.png"