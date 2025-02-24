#!/bin/sh

while getopts d:a:p:f:r:m:g: flag
do
    case "${flag}" in
        d) directory_new=${OPTARG};;  # assign -d to directory_new
        a) adapter_file=${OPTARG};;
        p) paired_end=${OPTARG};;
        f) pre_filtered=${OPTARG};;
        r) reference_directory=${OPTARG};;
        m) miRBase_species_abbreviation=${OPTARG};;
        g) genome_species_abbreviation=${OPTARG};;
    esac
done

if [ -z "$directory_new" ]; then
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

if [ -z "$pre_filtered" ]; then
    echo "Error: Filter out pre-miRNAs and tRNAs? [-f]"
    exit 1
fi

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

#Add / to the end of path to directory_new if user did not include one
directory_new=$(echo "${directory_new}" | sed 's![^/]$!&/!')

#Add / to the end of path to reference directory_new if user did not include one
reference_directory=$(echo "${reference_directory}" | sed 's![^/]$!&/!')

#Parameters for filtering pre-miRNA BLAST results
minimum_pre_miRNA_evalue=.05
minimum_pre_miRNA_bit_score=50

#Parameters for filtering tRNA BLAST results
minimum_tRNA_evalue=.05
minimum_tRNA_bit_score=50

#Parameters for filtering BLAST results
minimum_evalue=0.1 #.05
minimum_length=14
maximum_mismatch=2
maximum_gap=1
maximum_gap_if_maximum_mismatch=0

#Minimum length of sequence after the sncRNA within a read (after adapter and barcode removal)
minimum_length_after_sncRNA=15

#File is generated in sample folder called "Sample".summary.txt that summarizes number of reads following each step
#Sections separated with ########## are code for counting reads

##############################################################################
#Activate conda environment
	location=$(conda info | awk '/base environment/' | awk '{print $4}')
	source ${location}/etc/profile.d/conda.sh

	conda activate SCRAP

samples=$(awk '!/^#/' ${adapter_file} | awk '{print $1}')

###Quality Control (FastQC)###
#Make output directory_new for FastQC reports
	mkdir ${directory_new}FastQC_Reports

for sample in $samples
do

  #Run FastQC
  #For further usage details: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/
  for file in $(ls ${directory_new}${sample}/ | awk '/fastq/')
  do

    fastqc \
    -o ${directory_new}FastQC_Reports \
    ${directory_new}${sample}/${file}

  done

done

###Quality Control (MultiQC)###
#Make output directory_new for MultiQC report
	mkdir ${directory_new}MultiQC_Report

#Run MultiQC
#For further usage details: https://multiqc.info/docs/
	multiqc ${directory_new}FastQC_Reports -o ${directory_new}MultiQC_Report

for sample in $samples
do
#Extract adapter and barcode sequences from adapter file
five_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $2}')
three_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $3}')
five_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $4}')
three_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $5}')

#Remove 'sample.summary.txt' file if it exists
	[ -e ${directory_new}${sample}/${sample}.summary.txt ] && rm ${directory_new}${sample}/${sample}.summary.txt

#Print the name of the sample being processed in the 'sample.summary.txt' file
	echo "$sample" >> ${directory_new}${sample}/${sample}.summary.txt

#Print the structure of the reads for the sample being processed in the 'sample.summary.txt' file
	echo "5'-${five_prime_adapter}${five_prime_barcode}...sncRNA-targetRNA...${three_prime_barcode}${three_prime_adapter}-3'" >> ${directory_new}${sample}/${sample}.summary.txt

#Print the miRBase species abbreviation used to process the sample in the 'sample.summary.txt' file
	echo "miRBase Species Abbreviation: ${miRBase_species_abbreviation}" >> ${directory_new}${sample}/${sample}.summary.txt

#Print the genome species abbreviation used to process the sample in the 'sample.summary.txt' file
	echo "Genome Species Abbreviation: ${genome_species_abbreviation}" >> ${directory_new}${sample}/${sample}.summary.txt

#Print the date and time that processing begins on the sample being processed in the 'sample.summary.txt' file
	echo "Start: $(date)" >> ${directory_new}${sample}/${sample}.summary.txt

#List all of the files in the sample folder (1 file for single-end sequencing, 2 files for paired-end sequencing)
for file in $(ls ${directory_new}${sample}/ | awk '/fastq/')
do

##########
#Print the number of reads in each of the files in the sample directory_new in the 'sample.summary.txt' file
	gunzip -c ${directory_new}${sample}/${file} > ${directory_new}${sample}/${file}.tmp
	echo "${file} raw reads: $(( $(wc -l ${directory_new}${sample}/${file}.tmp | awk '{print $1}') / 4 ))" >> ${directory_new}${sample}/${sample}.summary.txt
	rm ${directory_new}${sample}/${file}.tmp
##########

done

###Combine Paired-End Reads (FLASh)###

#If the user set the paired-end parameter [-p] to yes, combine paired end reads
#IF the user set the paired-end parameter [-p] to no (or anything other than yes), skip this step
if [ $paired_end == "yes" ]
then

#Make directory_new for FLASh output
	mkdir ${directory_new}${sample}/${sample}_FLASh

#Run FLASh
#For further usage details: http://gensoft.pasteur.fr/docs/FLASH/1.2.11/flash
	flash \
	--allow-outies \
	--output-directory_new=${directory_new}${sample}/${sample}_FLASh/ \
	--output-prefix=${sample} \
	--max-overlap=150 \
	--min-overlap=6 \
	--compress \
	${directory_new}${sample}/${sample}_R1.fastq.gz \
	${directory_new}${sample}/${sample}_R2.fastq.gz \
	2>&1 | tee ${directory_new}${sample}/${sample}_FLASh/FLASh_${sample}.log

#Move and rename the output of FLASh to '.fastq.gz' so that input file name for Cutadapt is the same whether starting with single end or paired-end reads
	mv ${directory_new}${sample}/${sample}_FLASh/${sample}.extendedFrags.fastq.gz ${directory_new}${sample}/${sample}.fastq.gz

#If paired-end reads were combined, delete the folder containing the temporary files generated during FLASh
	[ -e ${directory_new}${sample}/${sample}_FLASh ] && rm -r ${directory_new}${sample}/${sample}_FLASh

##########
#Print the number of reads following FLASh (if performed) in the 'sample.summary.txt' file
	gunzip -k ${directory_new}${sample}/${sample}.fastq.gz
	echo "${sample} combined paired-end reads: $(( $(wc -l ${directory_new}${sample}/${sample}.fastq | awk '{print $1}') / 4 ))" >> ${directory_new}${sample}/${sample}.summary.txt
	rm ${directory_new}${sample}/${sample}.fastq
##########

fi

#Copy the output of FLASh to a new file so that the original file does not get deleted
	cp ${directory_new}${sample}/${sample}.fastq.gz ${directory_new}${sample}/${sample}.tmp.fastq.gz

###Remove Adapters without Random Nucleotides, Qualty Trimming, Minimum Length Threshold (Cutadapt)###

#Run Cutadapt
#For further usage details: https://cutadapt.readthedocs.io/en/stable/guide.html

#If present, remove the 5' adapter sequence from the read
#Minimum quality score set to 30 [-q 30]
#Minimum length set to 30 [-m 30]
#Some library preparations may result in sequential ligation of adapters, therefore 2 rounds [-n 2] of adapter removal are performed in case this occurred

if [ ! -z "${five_prime_adapter}" ]
then

	cutadapt \
	-g ${five_prime_adapter} \
	-q 30 \
	-m 30 \
	-n 2 \
	-j 12 \
	-o ${directory_new}${sample}/${sample}.cutadapt.fastq.gz \
	--json=${directory_new}${sample}/${sample}.cutadapt.5adapter.json \
	${directory_new}${sample}/${sample}.tmp.fastq.gz

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed
	mv ${directory_new}${sample}/${sample}.cutadapt.fastq.gz ${directory_new}${sample}/${sample}.tmp.fastq.gz

fi

#If present, remove the 3' adapter sequence from the read
#Parameters are the same as for the 5' adapter sequence removal
if [ ! -z "${three_prime_adapter}" ]
then

	cutadapt \
	-a ${three_prime_adapter} \
	-q 30 \
	-m 30 \
	-n 2 \
	-j 12 \
	-o ${directory_new}${sample}/${sample}.cutadapt.fastq.gz \
	--json=${directory_new}${sample}/${sample}.cutadapt.3adapter.json \
	${directory_new}${sample}/${sample}.tmp.fastq.gz

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed
	mv ${directory_new}${sample}/${sample}.cutadapt.fastq.gz ${directory_new}${sample}/${sample}.tmp.fastq.gz

fi

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed
	mv ${directory_new}${sample}/${sample}.tmp.fastq.gz ${directory_new}${sample}/${sample}.cutadapt.fastq.gz

# Create a directory for FastQC post-trimming results if it doesn't exist
mkdir -p ${directory_new}${sample}/fastQC_post_trimming

# Run FastQC on the new Cutadapt output and store results in the new directory
fastqc -t 12 -o ${directory_new}${sample}/fastQC_post_trimming/ ${directory_new}${sample}/${sample}.cutadapt.fastq.gz

# Trim reads to a maximum length of 80bp using Cutadapt (use directory_new consistently)
cutadapt -l 80 -o ${directory_new}${sample}/${sample}.trimmed.fastq.gz ${directory_new}${sample}/${sample}.cutadapt.fastq.gz

# Ensure output maintains correct name for downstream steps
mv ${directory_new}${sample}/${sample}.trimmed.fastq.gz ${directory_new}${sample}/${sample}.cutadapt.fastq.gz

###Deduplicate Reads###

#Unzip FASTQ file from Cutadapt
	gunzip ${directory_new}${sample}/${sample}.cutadapt.fastq.gz

##########
#Print the number of reads following nonbarcoded adapter removal in the 'sample.summary.txt' file
	echo "${sample} reads following unbarcoded adapter removal: $(( $(wc -l ${directory_new}${sample}/${sample}.cutadapt.fastq | awk '{print $1}') / 4 ))" >> ${directory_new}${sample}/${sample}.summary.txt
##########

#Remove duplicate reads, and changes the read ID
#Read ID now contains rank (by duplicate number) and number of duplicates
	awk 'NR%4==2' ${directory_new}${sample}/${sample}.cutadapt.fastq | \
	sort -T ${directory_new}${sample} | \
	uniq -c | \
	sort -k1,1nr -T ${directory_new}${sample} | \
	awk '{print $0,NR}' | \
	awk '{print ">"$3"-"$1"\n"$2}' \
	> ${directory_new}${sample}/${sample}.cutadapt.deduped.fasta

	rm ${directory_new}${sample}/${sample}.cutadapt.fastq

##########
#Print the number of reads following deduplication in the 'sample.summary.txt' file
	echo "${sample} deduplicated reads: $(( $(wc -l ${directory_new}${sample}/${sample}.cutadapt.deduped.fasta | awk '{print $1}') / 2 ))" >> ${directory_new}${sample}/${sample}.summary.txt
##########

###Remove random barcodes (Cutadapt)###

#If present, remove the 5' barcode sequence from the read
#Parameters are the same as for the 5' adapter sequence removal
if [ ! -z "${five_prime_barcode}" ]
then

	cutadapt \
	-g ^${five_prime_barcode} \
	-m 30 \
	-j 12 \
	-o ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	--json=${directory_new}${sample}/${sample}.cutadapt.5barcode.json \
	${directory_new}${sample}/${sample}.cutadapt.deduped.fasta

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed
	mv ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta ${directory_new}${sample}/${sample}.cutadapt.deduped.fasta
fi

#If present, remove the 3' barcode sequence from the read
#Parameters are the same as for the 5' adapter sequence removal
if [ ! -z "${three_prime_barcode}" ]
then

	cutadapt \
	-a ${three_prime_barcode}$ \
	-m 30 \
	-j 12 \
	-o ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	--json=${directory_new}${sample}/${sample}.cutadapt.3barcode.json \
	${directory_new}${sample}/${sample}.cutadapt.deduped.fasta

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed
	mv ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta ${directory_new}${sample}/${sample}.cutadapt.deduped.fasta
fi

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed
	mv ${directory_new}${sample}/${sample}.cutadapt.deduped.fasta ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta

##########

#Print the number of reads following barcoded adapter removal in the 'sample.summary.txt' file
	echo "${sample} reads following barcode removal: $(( $(wc -l ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta | awk '{print $1}') / 2 ))" >> ${directory_new}${sample}/${sample}.summary.txt
##########

###Filter our pre-miRNAs and tRNAs (BLAST)###
if [ $pre_filtered == "yes" ]
then

#Align reads to pre-miRNA reference file
	blastn \
	-db ${reference_directory}/fasta/${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
	-query ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	-out ${directory_new}${sample}/${sample}.preprocessed.premiRNA.blast \
	-word_size 11 \
	-outfmt 6 \
	-num_threads 1 \
	-strand plus

#Filter pre-miRNA alignments
	awk -v var1=$minimum_pre_miRNA_evalue '($11 + 0) <= var1' ${directory_new}${sample}/${sample}.preprocessed.premiRNA.blast | \
	awk -v var1=$minimum_pre_miRNA_bit_score '$12 >= var1' | \
	awk '!x[$1]++' \
	> ${directory_new}${sample}/${sample}.preprocessed.premiRNA.blast.filtered

	rm ${directory_new}${sample}/${sample}.preprocessed.premiRNA.blast

#Remove reads that aligned to pre-miRNAs with minimum filters
	join -j1 -v2 \
	<(awk '{print $1}' ${directory_new}${sample}/${sample}.preprocessed.premiRNA.blast.filtered | sort -k1,1 -T ${directory_new}${sample} | uniq) \
	<(awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta | sort -k1,1 -T ${directory_new}${sample}) | \
	awk '{print ">"$1"\n"$2}' \
	> ${directory_new}${sample}/${sample}.premiRNA.removed.fasta

##########
#Print the number of reads following pre-miRNA filtering in the 'sample.summary.txt' file
	echo "${sample} pre-miRNA filtered: $(( $(wc -l ${directory_new}${sample}/${sample}.premiRNA.removed.fasta | awk '{print $1}') / 2 ))" >> ${directory_new}${sample}/${sample}.summary.txt
##########

	rm ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta
	rm ${directory_new}${sample}/${sample}.preprocessed.premiRNA.blast.filtered

#Align reads to tRNA reference file
	blastn \
	-db ${reference_directory}/fasta/${miRBase_species_abbreviation}/tRNA_${miRBase_species_abbreviation}.fasta \
	-query ${directory_new}${sample}/${sample}.premiRNA.removed.fasta \
	-out ${directory_new}${sample}/${sample}.preprocessed.tRNA.blast \
	-word_size 11 \
	-outfmt 6 \
	-num_threads 1 \
	-strand plus

#Filter tRNA alignments
	awk -v var1=$minimum_tRNA_evalue '($11 + 0) <= var1' ${directory_new}${sample}/${sample}.preprocessed.tRNA.blast | \
	awk -v var1=$minimum_tRNA_bit_score '$12 >= var1' | \
	awk '!x[$1]++' \
	> ${directory_new}${sample}/${sample}.preprocessed.tRNA.blast.filtered

	rm ${directory_new}${sample}/${sample}.preprocessed.tRNA.blast

#Remove reads that aligned to tRNAs with minimum filters
	join -j1 -v2 \
	<(awk '{print $1}' ${directory_new}${sample}/${sample}.preprocessed.tRNA.blast.filtered | sort -k1,1 -T ${directory_new}${sample} | uniq) \
	<(awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' ${directory_new}${sample}/${sample}.premiRNA.removed.fasta | sort -k1,1 -T ${directory_new}${sample}) | \
	awk '{print ">"$1"\n"$2}' \
	> ${directory_new}${sample}/${sample}.tRNA.removed.fasta

##########
#Print the number of reads following tRNA filtering in the 'sample.summary.txt' file
	echo "${sample} tRNA filtered: $(( $(wc -l ${directory_new}${sample}/${sample}.tRNA.removed.fasta | awk '{print $1}') / 2 ))" >> ${directory_new}${sample}/${sample}.summary.txt
##########

	rm ${directory_new}${sample}/${sample}.premiRNA.removed.fasta
	rm ${directory_new}${sample}/${sample}.preprocessed.tRNA.blast.filtered

#Rename the output of filtering so that it matches whether or not filtering was performed
	mv ${directory_new}${sample}/${sample}.tRNA.removed.fasta ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta

fi

###Identify sncRNA (BLAST)###

#Run BLAST
#For further usage details: https://www.ncbi.nlm.nih.gov/books/NBK279690/pdf/Bookshelf_NBK279690.pdf

#Align reads to sncRNA reference (configured during installation)
	makeblastdb \
	-in ${reference_directory}/fasta/${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
	-dbtype nucl

	blastn \
	-db ${reference_directory}/fasta/${miRBase_species_abbreviation}/sncRNA_${miRBase_species_abbreviation}.fasta \
	-query ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	-out ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
	-word_size 11 \
	-outfmt 6 \
	-num_threads 12

	mv ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.sncRNA.blast

	# BLAST the filtered CLIPPED reads with the miRNA sequences (plus strand)
	blastn \
	-db ${reference_directory}/fasta/${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
	-query ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	-out ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
	-word_size 11 \
	-outfmt 6 \
	-num_threads 12 \
	-strand plus

	mv ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.plus.blast

	# BLAST the filtered CLIPPED reads with the miRNA sequences (minus strand)
	blastn \
	-db ${reference_directory}/fasta/${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
	-query ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	-out ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
	-word_size 11 \
	-outfmt 6 \
	-num_threads 12 \
	-strand minus

	mv ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.minus.blast

	# BLAST the miRNA with the hairpin sequences
	blastn \
	-db ${reference_directory}/fasta/${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
	-query ${reference_directory}/fasta/${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
	-out ${directory_new}${sample}/${sample}.hairpin.miRNA.blast \
	-word_size 11 \
	-outfmt 6 \
	-num_threads 12

	# BLAST the filtered CLIPPED reads with the hairpin sequences (plus strand)
	blastn \
	-db ${reference_directory}/fasta/${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
	-query ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	-out ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
	-word_size 11 \
	-outfmt 6 \
	-num_threads 12 \
	-strand plus

	mv ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.plus.blast

	# BLAST the filtered CLIPPED reads with the hairpin sequences (minus strand)
	blastn \
	-db ${reference_directory}/fasta/${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
	-query ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	-out ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
	-word_size 11 \
	-outfmt 6 \
	-num_threads 12 \
	-strand minus

	mv ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.minus.blast

#Filter BLAST results (parameters indicated at the beginning of the script)
	awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.plus.blast | \
	awk -v var1=$minimum_length '$4 >= var1' | \
	awk -v var1=$maximum_mismatch '$5 <= var1' | \
	awk -v var1=$maximum_gap '$6 <= var1' | \
	awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
	awk '!x[$1]++' \
	> ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.plus.blast.filtered

	awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.minus.blast | \
	awk -v var1=$minimum_length '$4 >= var1' | \
	awk -v var1=$maximum_mismatch '$5 <= var1' | \
	awk -v var1=$maximum_gap '$6 <= var1' | \
	awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
	awk '!x[$1]++' \
	> ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.minus.blast.filtered

	awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.plus.blast | \
	awk -v var1=$minimum_length '$4 >= var1' | \
	awk -v var1=$maximum_mismatch '$5 <= var1' | \
	awk -v var1=$maximum_gap '$6 <= var1' | \
	awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
	awk '!x[$1]++' \
	> ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.plus.blast.filtered

	awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.minus.blast | \
	awk -v var1=$minimum_length '$4 >= var1' | \
	awk -v var1=$maximum_mismatch '$5 <= var1' | \
	awk -v var1=$maximum_gap '$6 <= var1' | \
	awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
	awk '!x[$1]++' \
	> ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.minus.blast.filtered

	join \
	<(sort -k1,1 ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.plus.blast.filtered -T ${directory_new}${sample}) \
	<(sort -k1,1 ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.plus.blast.filtered -T ${directory_new}${sample}) \
	> ${directory_new}${sample}/${sample}.blast.miRNA_iso.plus.merged

	join \
	<(sort -k1,1 ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.minus.blast.filtered -T ${directory_new}${sample}) \
	<(sort -k1,1 ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.minus.blast.filtered -T ${directory_new}${sample}) \
	> ${directory_new}${sample}/${sample}.blast.miRNA_iso.minus.merged

	# joined on the file with the miRNA alignments on the hairpin
	join -1 13 -2 2 \
	<(sort -k13 ${directory_new}${sample}/${sample}.blast.miRNA_iso.plus.merged) \
	<(sort -k2 ${directory_new}${sample}/${sample}.hairpin.miRNA.blast) \
	| awk '{if($3==$24)print $0}' | sort -k2 > temp1.plus

	join -1 13 -2 2 \
	<(sort -k13 ${directory_new}${sample}/${sample}.blast.miRNA_iso.minus.merged) \
	<(sort -k2 ${directory_new}${sample}/${sample}.hairpin.miRNA.blast) \
	| awk '{if($3==$24)print $0}' | sort -k2 > temp1.minus

	####################################################################################

# The objective is to get the isomiRic miRNA alignments
# joined on the file with the miRNA alignments on the hairpin
	join -1 13 -2 2 \
	<(sort -k13 ${directory_new}${sample}/${sample}.blast.miRNA_iso.plus.merged) \
	<(sort -k2 ${directory_new}${sample}/${sample}.hairpin.miRNA.blast) \
	| awk '{if($3==$24)print $0}' | sort -k2 > temp1.plus

	join -1 13 -2 2 \
	<(sort -k13 ${directory_new}${sample}/${sample}.blast.miRNA_iso.minus.merged) \
	<(sort -k2 ${directory_new}${sample}/${sample}.hairpin.miRNA.blast) \
	| awk '{if($3==$24)print $0}' | sort -k2 > temp1.minus

	# check if the alignment length is different:
awk 'function abs(x){return ((x < 0.0) ? -x : x)} \
		{ \
			if($26 >= 17){
			# check on the direction of the miRNA aligning to hairpin
				if($32 > $31){
					if($32 > $21 && $31 > $20) \
						print $0,"\t3p:-"abs($32-$21)";5p:+"abs($31-$20);
					else if($32 == $21 && $31 > $20) \
						print $0,"\t5p:+"abs($31-$20);
					else if($32 == $21 && $31 < $20) \
						print $0,"\t5p:-"abs($31-$20);
					else if($32 < $21 && $31 > $20) \
						print $0 ,"\t3p:+"abs($32-$21)";5p:+"abs($31-$20);
					else if($32 < $21 && $31 == $20) \
						print $0 ,"\t3p:+"abs($32-$21);
					else if($32 > $21 && $31 < $20) \
						print $0,"\t3p:-"abs($32-$21)";5p:-"abs($31-$20);
					else if($32 > $21 && $31 == $20) \
						print $0,"\t3p:-"abs($32-$21);
					else if($32 < $21 && $31 < $20) \
						print $0 ,"\t3p:+"abs($32-$21)";5p:-"abs($31-$20);
					else \
						print $0; \
				}
				# else
					# print $0"\tND: miRNA alignment direction";
			}
			else
				print $0"\tND\tmiRNA alignment length";
		}' temp1.plus > plus #| less # > SRR959751.blast.possible.miRNA_iso

awk 'function abs(x){return ((x < 0.0) ? -x : x)} \
		{ \
			if($26 >= 17){
			# check on the direction of the miRNA aligning to hairpin
				if($32 > $31){
					if($32 > $20 && $31 > $21) \
						print $0,"\t3p:-"abs($32-$20)";5p:+"abs($31-$21);

					else if($32 == $20 && $31 > $21) \
						print $0,"\t5p:+"abs($31-$21);

					else if($32 == $20 && $31 < $21) \
						print $0,"\t5p:-"abs($31-$21);

					else if($32 < $20 && $31 > $21) \
						print $0 ,"\t3p:+"abs($32-$20)";5p:+"abs($31-$21);

					else if($32 < $20 && $31 == $21) \
						print $0 ,"\t3p:+"abs($32-$20);

					else if($32 > $20 && $31 < $21) \
						print $0,"\t3p:-"abs($32-$20)";5p:-"abs($31-$21);

					else if($32 > $20 && $31 == $21) \
						print $0,"\t3p:-"abs($32-$20);

					else if($32 < $20 && $31 < $21) \
						print $0 ,"\t3p:+"abs($32-$20)";5p:-"abs($31-$21);

					else \
						print $0; \
				}
				# else
				# 	print $0"\tND\tmiRNA alignment direction";
			}
			else
				print $0"\tND\tmiRNA alignment length";
		}' temp1.minus > minus #| less # > SRR959751.blast.possible.miRNA_iso

	# pull out the canonical miRNA alignments from both files, and make the alignment file for further processing
	awk '{if(NF==34)print $2"\t"$3"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' plus > plus.miRNA
	awk '{if(NF==34)print $2"\t"$3"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' minus > minus.miRNA

	# pull out the isomiRic miRNA alignments from both files
	awk '{if(NF==35)print $2"\t"$3"|"$NF"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' plus > plus.miRNA_isomiR
	awk '{if(NF==35)print $2"\t"$3"|"$NF"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' minus > minus.miRNA_isomiR

	# put the resulting alignment in the total annotation file
	cat plus.miRNA minus.miRNA plus.miRNA_isomiR minus.miRNA_isomiR > test.miRNA.alignments.blast

####################################################################################
#Make tabular version of barcoded FASTA file (from before BLAST)
	awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' \
	${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	> ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.tab

# rename the file to join the pipeline
mv test.miRNA.alignments.blast ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered

#=====================================================================================================
#Join tabular input file with filtered BLAST results
	join \
	<(sort -k1,1 ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered -T ${directory_new}${sample}) \
	<(sort -k1,1 ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.tab -T ${directory_new}${sample}) \
	> ${directory_new}${sample}/${sample}.blast.merged

	rm ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered
	rm ${directory_new}${sample}/${sample}.cutadapt.deduped.barcoded.tab

##########
#Print the number of reads containing a sncRNA (satisfying BLAST filters) in the 'sample.summary.txt' file
	echo "${sample} filtered BLAST alignments: $(( $(wc -l ${directory_new}${sample}/${sample}.blast.merged | awk '{print $1}') / 1 ))" >> ${directory_new}${sample}/${sample}.summary.txt
##########
	echo "Ran successfully till here !!"
###Identify potential chimeras###

#Create FASTA file with sequence of sncRNA remove and name appended to the read ID
#Only output reads containing minimum length of sequence (parameter indicated at the beginning of the script) following sncRNA removal
	awk '{stop=$8;readlength=length($13);print $0,"\t", readlength-(stop)}' ${directory_new}${sample}/${sample}.blast.merged | \
	awk -v var1=$minimum_length_after_sncRNA '$14 >= var1' | \
	awk '{print $0,"\t",substr($13, $8+1, length($13)-($8))}' | \
	awk '{print">"$1"."$2"\n"$15}' \
	> ${directory_new}${sample}/${sample}.target.fasta

	rm ${directory_new}${sample}/${sample}.blast.merged

##########
#Print the number of reads containing a sncRNA (satisfying BLAST filters) and minimum length of sequence in the 'sample.summary.txt' file
	echo "${sample} potential chimeras: $(( $(wc -l ${directory_new}${sample}/${sample}.target.fasta | awk '{print $1}') / 2 ))" >> ${directory_new}${sample}/${sample}.summary.txt
##########

###Align reads to the reference genome (HISAT2)##
#For further usage details: http://daehwankimlab.github.io/hisat2/manual/

#Reference configured during installation
#Run HISAT2
hisat2 \
	-x ${reference_directory}/fasta/${genome_species_abbreviation}/${genome_species_abbreviation} \
	-f ${directory_new}${sample}/${sample}.target.fasta \
	-S ${directory_new}${sample}/${sample}.aligned.sam \
	--summary-file ${directory_new}${sample}/${sample}.hisat2summary.txt
##########
#Print the number of reads that aligned the genome a single time in the 'sample.summary.txt' file

	echo "${sample} unique alignments: $(( $(awk '$5 == 60' ${directory_new}${sample}/${sample}.aligned.sam | wc -l | awk '{print $1}') / 1 ))" >> ${directory_new}${sample}/${sample}.summary.txt
##########

#Select only uniquely mapped reads

	awk '/^@/ || $5 == 60' ${directory_new}${sample}/${sample}.aligned.sam \
	> ${directory_new}${sample}/${sample}.aligned.unique.sam

#Convert to bam file

	samtools view \
	-S \
	-h \
	-b ${directory_new}${sample}/${sample}.aligned.unique.sam \
	> ${directory_new}${sample}/${sample}.aligned.unique.bam

echo "Finish: $(date)" >> ${directory_new}${sample}/${sample}.summary.txt

done

conda deactivate
