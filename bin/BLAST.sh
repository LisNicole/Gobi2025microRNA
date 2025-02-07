#!/bin/bash

# Parse command-line arguments
while getopts r:m:g:s: flag
do
    case "${flag}" in
        r) reference_directory=${OPTARG};;
        m) miRBase_species_abbreviation=${OPTARG};;
        g) genome_species_abbreviation=${OPTARG};;
        s) species=${OPTARG};;
    esac
done

# Check if required arguments are provided
if [ -z "$reference_directory" ]; then
    echo "Error: Path to reference directory not provided [-r]"
    exit 1
fi

if [ -z "$miRBase_species_abbreviation" ]; then
    echo "Error: miRBase species abbreviation not provided [-m]"
    exit 1
fi

if [ -z "$genome_species_abbreviation" ]; then
    echo "Error: Species genome abbreviation not provided [-g]"
    exit 1
fi

if [ -z "$species" ]; then
    echo "Error: Species not provided [-s]"
    exit 1
fi

# Change to the reference directory
cd ${reference_directory} || exit
cd fasta || exit

# Create directories for miRBase and genome species
mkdir -p ${miRBase_species_abbreviation}
mkdir -p ${genome_species_abbreviation}

# Subset miRNA sequences for species of interest
# miRBase.fasta obtained using: wget https://www.mirbase.org/download/mature.fa -O miRBase.fasta
# tRFdb.fasta obtained using: github from SCRAP
# Final miRNA name format (example): miRNA-mmu-let-7c-5p

awk '/^>/ {printf "%s%s ", pfx, $0; pfx="\n"; next} {printf "%s", $0} END {print ""}' miRBase.fasta | \
grep "${miRBase_species_abbreviation}-" | \
awk '{print $1"\n"$6}' | \
sed 's/>/>miRNA-/g' \
> ${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta

grep -A1 "${miRBase_species_abbreviation}" tRFdb.fasta > ${miRBase_species_abbreviation}/tRF_${miRBase_species_abbreviation}.fasta

cat ${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
    ${miRBase_species_abbreviation}/tRF_${miRBase_species_abbreviation}.fasta \
> ${miRBase_species_abbreviation}/sncRNA_${miRBase_species_abbreviation}.fasta

# Make BLAST database for sncRNA list
makeblastdb -in ${miRBase_species_abbreviation}/sncRNA_${miRBase_species_abbreviation}.fasta -dbtype nucl

# Download genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz -O mm39.fa.gz
gunzip mm39.fa.gz

# Index mouse genome for HISAT2 alignment (takes ~2 hours)
cd ${genome_species_abbreviation} || exit

hisat2-build /u/halle/ge59puj/home_at/gobi2025/reference/fasta/mm39.fa mm39

cd ..
# Download reference genome chromosome sizes
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm39/bigZips/mm39.chrom.sizes -O mm39.chrom.sizes

#Download refFlat

	wget https://hgdownload.cse.ucsc.edu/goldenPath/mm39/database/refFlat.txt.gz -O mm39refFlat.txt

#Index mouse genome for GATK and Samtools


	gatk CreateSequenceDictionary -R mm39/mm39.fa
	samtools faidx mm39/mm39.fa

#Subset miRNA hairpin sequences for species of interest
#miRBase.hairpin.fasta obtained using:	wget https://www.mirbase.org/ftp/CURRENT/hairpin.fa.gz -O miRBase.hairpin.fasta

	awk '/^>/ {printf "%s%s ", pfx, $0; pfx="\n"; next} {printf "%s", $0} END {print ""}' miRBase.hairpin.fasta | \
	grep "mm39-" | \
	awk '{print $1"\n"$7}' | \
	sed 's/>/>miRNA-hairpin-/g' \
	> mm39/hairpin_mm39.fasta

	#Make BLAST database for miRNA hairpin sequences

	makeblastdb \
	-in mm39/hairpin_mm39.fasta \
	-dbtype nucl


#Subset tRNA sequences for species of interest

	grep -A1 ">mm39_" GtRNAdb.fasta \
	> tRNA_mm39.fasta

#Make BLAST database for miRNA hairpin sequences

	makeblastdb \
	-in tRNA_mm39.fasta \
	-dbtype nucl

	cd ..
	cd annotation

  gunzip ${species}.annotation.bed.gz

	conda deactivate
