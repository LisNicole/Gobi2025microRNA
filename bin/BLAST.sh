#!/bin/sh

# Parse arguments
while getopts r:m:g:s: flag
do
    case "${flag}" in
        r) reference_directory=${OPTARG};;
        m) miRBase_species_abbreviation=${OPTARG};;
        g) genome_species_abbreviation=${OPTARG};;
        s) species=${OPTARG};;
    esac
done

# Check if necessary arguments are provided
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

# Check if virtual environment is activated
if [ -z "$VIRTUAL_ENV" ]; then
    echo "Error: Virtual environment is not activated. Please activate your virtual environment first."
    exit 1
fi

# Go to the reference directory and create necessary directories
cd ${reference_directory} || { echo "Reference directory not found!"; exit 1; }

mkdir -p ${miRBase_species_abbreviation}
mkdir -p ${genome_species_abbreviation}
mkdir -p annotation  # Create annotation directory if it doesn't exist

# Download the missing reference files if not already present
if [ ! -f "miRBase.fasta" ]; then
    echo "Downloading miRBase.fasta..."
    wget -O miRBase.fasta http://example.com/path/to/miRBase.fasta  # Replace with actual URL
fi

if [ ! -f "tRFdb.fasta" ]; then
    echo "Downloading tRFdb.fasta..."
    wget -O tRFdb.fasta http://example.com/path/to/tRFdb.fasta  # Replace with actual URL
fi

if [ ! -f "miRBase.hairpin.fasta" ]; then
    echo "Downloading miRBase.hairpin.fasta..."
    wget -O miRBase.hairpin.fasta http://example.com/path/to/miRBase.hairpin.fasta  # Replace with actual URL
fi

if [ ! -f "GtRNAdb.fasta" ]; then
    echo "Downloading GtRNAdb.fasta..."
    wget -O GtRNAdb.fasta http://example.com/path/to/GtRNAdb.fasta  # Replace with actual URL
fi

# Create BLAST databases if necessary files exist
if [ -f "miRBase.fasta" ]; then
    awk '/^>/ {printf "%s%s ", pfx, $0; pfx="\n"; next} {printf "%s", $0} END {print ""}' miRBase.fasta | \
    grep "${miRBase_species_abbreviation}-" | \
    awk '{print $1"\n"$6}' | \
    sed 's/>/>miRNA-/g' \
    > ${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta

    grep -A1 "${miRBase_species_abbreviation}" tRFdb.fasta \
    > ${miRBase_species_abbreviation}/tRF_${miRBase_species_abbreviation}.fasta

    cat \
    ${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
    ${miRBase_species_abbreviation}/tRF_${miRBase_species_abbreviation}.fasta \
    > ${miRBase_species_abbreviation}/sncRNA_${miRBase_species_abbreviation}.fasta

    # Make BLAST database for sncRNA list
    makeblastdb \
    -in ${miRBase_species_abbreviation}/sncRNA_${miRBase_species_abbreviation}.fasta \
    -dbtype nucl
fi

# Subset miRNA hairpin sequences for species of interest
if [ -f "miRBase.hairpin.fasta" ]; then
    awk '/^>/ {printf "%s%s ", pfx, $0; pfx="\n"; next} {printf "%s", $0} END {print ""}' miRBase.hairpin.fasta | \
    grep "${miRBase_species_abbreviation}-" | \
    awk '{print $1"\n"$7}' | \
    sed 's/>/>miRNA-hairpin-/g' \
    > ${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta

    # Make BLAST database for miRNA hairpin sequences
    makeblastdb \
    -in ${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
    -dbtype nucl
fi

# Subset tRNA sequences for species of interest
if [ -f "GtRNAdb.fasta" ]; then
    grep -A1 ">${miRBase_species_abbreviation}_" GtRNAdb.fasta \
    > ${miRBase_species_abbreviation}/tRNA_${miRBase_species_abbreviation}.fasta

    # Make BLAST database for tRNA sequences
    makeblastdb \
    -in ${miRBase_species_abbreviation}/tRNA_${miRBase_species_abbreviation}.fasta \
    -dbtype nucl
fi

# Download and process the genome sequences for the specified species
echo "Downloading genome for ${genome_species_abbreviation}..."

mkdir -p ${genome_species_abbreviation}

# Download the reference genome sequence in fasta format if not already present
if [ ! -f "${genome_species_abbreviation}/${genome_species_abbreviation}.fa" ]; then
    wget http://hgdownload.cse.ucsc.edu/goldenpath/${genome_species_abbreviation}/bigZips/${genome_species_abbreviation}.fa.gz -O ${genome_species_abbreviation}/${genome_species_abbreviation}.fa.gz
    gunzip -c ${genome_species_abbreviation}/${genome_species_abbreviation}.fa.gz > ${genome_species_abbreviation}/${genome_species_abbreviation}.fa
    rm -r ${genome_species_abbreviation}/${genome_species_abbreviation}.fa.gz
fi

# Index genome with HISAT2 (Takes time)
echo "Indexing genome with HISAT2..."
if ! command -v hisat2-build &> /dev/null; then
    echo "Error: hisat2-build command not found. Please install HISAT2."
    exit 1
fi

cd ${genome_species_abbreviation}
hisat2-build ${genome_species_abbreviation}.fa ${genome_species_abbreviation}
cd ..

# Download the chromosome sizes for the reference genome if not already present
if [ ! -f "${genome_species_abbreviation}/${genome_species_abbreviation}.chrom.sizes" ]; then
    echo "Downloading chromosome sizes for ${genome_species_abbreviation}..."
    wget http://hgdownload.cse.ucsc.edu/goldenpath/${genome_species_abbreviation}/bigZips/${genome_species_abbreviation}.chrom.sizes -O ${genome_species_abbreviation}/${genome_species_abbreviation}.chrom.sizes
fi

# Now go to annotation directory and unzip the annotation file
cd annotation
if [ ! -f "${species}.annotation.bed.gz" ]; then
    echo "Error: ${species}.annotation.bed.gz file not found!"
    exit 1
fi
gunzip ${species}.annotation.bed.gz

echo "Processing complete!"
