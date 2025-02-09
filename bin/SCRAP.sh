#!/bin/sh

#Command to run this script: bash gobi2025/bin/SCRAP.sh \
                                # -d /u/halle/ge59puj/home_at/gobi2025/directory/ \
                                # -a /u/halle/ge59puj/home_at/gobi2025/adapters/test.txt \
                                # -p no \
                                # -f yes \
                                # -r gobi2025/ \
                                # -m mmu \
                                # -g mm39

while getopts d:a:p:f:r:m:g: flag
do
    case "${flag}" in
        d) directory=${OPTARG};;
        a) adapter_file=${OPTARG};;
        p) paired_end=${OPTARG};;
        f) pre_filtered=${OPTARG};;
        r) reference_directory=${OPTARG};;
        m) miRBase_species_abbreviation=${OPTARG};;
        g) genome_species_abbreviation=${OPTARG};;
    esac
done

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

#Add / to the end of path to directory if user did not include one

directory=$(echo "${directory}" | sed 's![^/]$!&/!')

#Add / to the end of path to reference directory if user did not include one

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

#In this pipeline, intermediate files are removed to save space on the computer's memory
#If you wish to keep intermediate files (and have sufficient storage on your computer), consider removing or hiding the lines of code beginning with 'rm'

#File is generated in sample folder called "Sample".summary.txt that summarizes number of reads following each step
#Sections separated with ########## are code for counting reads

#For the most part, code modifications are not required below this line
##############################################################################

#Activate conda environment

	location=$(conda info | awk '/base environment/' | awk '{print $4}')
	source ${location}/etc/profile.d/conda.sh

	conda activate SCRAP

samples=$(awk '!/^#/' ${adapter_file} | awk '{print $1}')


###Quality Control (FastQC)###

#Make output directory for FastQC reports

	mkdir ${directory}FastQC_Reports

for sample in $samples
do

#Run FastQC
#For further usage details: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/

for file in $(ls ${directory}${sample}/ | awk '/fastq/')
do

	fastqc \
	-o ${directory}FastQC_Reports \
	${directory}${sample}/${file}

done

done


###Quality Control (MultiQC)###

#Make output directory for MultiQC report

	mkdir ${directory}MultiQC_Report

#Run MultiQC
#For further usage details: https://multiqc.info/docs/

	multiqc ${directory}FastQC_Reports -o ${directory}MultiQC_Report


for sample in $samples
do

#Extract adapter and barcode sequences from adapter file

five_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $2}')
three_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $3}')
five_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $4}')
three_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $5}')

#Remove 'sample.summary.txt' file if it exists

	[ -e ${directory}${sample}/${sample}.summary.txt ] && rm ${directory}${sample}/${sample}.summary.txt

#Print the name of the sample being processed in the 'sample.summary.txt' file

	echo "$sample" >> ${directory}${sample}/${sample}.summary.txt

#Print the structure of the reads for the sample being processed in the 'sample.summary.txt' file

	echo "5'-${five_prime_adapter}${five_prime_barcode}...sncRNA-targetRNA...${three_prime_barcode}${three_prime_adapter}-3'" >> ${directory}${sample}/${sample}.summary.txt

#Print the miRBase species abbreviation used to process the sample in the 'sample.summary.txt' file

	echo "miRBase Species Abbreviation: ${miRBase_species_abbreviation}" >> ${directory}${sample}/${sample}.summary.txt

#Print the genome species abbreviation used to process the sample in the 'sample.summary.txt' file

	echo "Genome Species Abbreviation: ${genome_species_abbreviation}" >> ${directory}${sample}/${sample}.summary.txt

#Print the date and time that processing begins on the sample being processed in the 'sample.summary.txt' file

	echo "Start: $(date)" >> ${directory}${sample}/${sample}.summary.txt

#List all of the files in the sample folder (1 file for single-end sequencing, 2 files for paired-end sequencing)

for file in $(ls ${directory}${sample}/ | awk '/fastq/')
do

##########
#Print the number of reads in each of the files in the sample directory in the 'sample.summary.txt' file

	gunzip -c ${directory}${sample}/${file} > ${directory}${sample}/${file}.tmp
	echo "${file} raw reads: $(( $(wc -l ${directory}${sample}/${file}.tmp | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
	rm ${directory}${sample}/${file}.tmp
##########

done


###Combine Paired-End Reads (FLASh)###

#If the user set the paired-end parameter [-p] to yes, combine paired end reads
#IF the user set the paired-end parameter [-p] to no (or anything other than yes), skip this step

if [ $paired_end == "yes" ]
then

#Make directory for FLASh output

	mkdir ${directory}${sample}/${sample}_FLASh

#Run FLASh
#For further usage details: http://gensoft.pasteur.fr/docs/FLASH/1.2.11/flash

	flash \
	--allow-outies \
	--output-directory=${directory}${sample}/${sample}_FLASh/ \
	--output-prefix=${sample} \
	--max-overlap=150 \
	--min-overlap=6 \
	--compress \
	${directory}${sample}/${sample}_R1.fastq.gz \
	${directory}${sample}/${sample}_R2.fastq.gz \
	2>&1 | tee ${directory}${sample}/${sample}_FLASh/FLASh_${sample}.log

#Move and rename the output of FLASh to '.fastq.gz' so that input file name for Cutadapt is the same whether starting with single end or paired-end reads

	mv ${directory}${sample}/${sample}_FLASh/${sample}.extendedFrags.fastq.gz ${directory}${sample}/${sample}.fastq.gz

#If paired-end reads were combined, delete the folder containing the temporary files generated during FLASh

	[ -e ${directory}${sample}/${sample}_FLASh ] && rm -r ${directory}${sample}/${sample}_FLASh


##########
#Print the number of reads following FLASh (if performed) in the 'sample.summary.txt' file

	gunzip -k ${directory}${sample}/${sample}.fastq.gz
	echo "${sample} combined paired-end reads: $(( $(wc -l ${directory}${sample}/${sample}.fastq | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
	rm ${directory}${sample}/${sample}.fastq
##########

fi

#Copy the output of FLASh to a new file so that the original file does not get deleted

	cp ${directory}${sample}/${sample}.fastq.gz ${directory}${sample}/${sample}.tmp.fastq.gz

###Remove Adapters without Random Nucleotides, Qualty Trimming, Minimum Length Threshold (Cutadapt)###

#Run Cutadapt
#For further usage details: https://cutadapt.readthedocs.io/en/stable/guide.html

#If present, remove the 5' adapter sequence from the read
#Minimum quality score set to 30 [-q 30]
#Minimum length set to 30 [-m 30]
#Some library preparations may result in sequential ligation of adapters, therefore 2 rounds [-n 2] of adapter removal are performed in case this occured

if [ ! -z "${five_prime_adapter}" ]
then

	cutadapt \
	-g ${five_prime_adapter} \
	-q 30 \
	-m 30 \
	-n 2 \
	-j 12 \
	-o ${directory}${sample}/${sample}.cutadapt.fastq.gz \
	--json=${directory}${sample}/${sample}.cutadapt.5adapter.json \
	${directory}${sample}/${sample}.tmp.fastq.gz

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed

	mv ${directory}${sample}/${sample}.cutadapt.fastq.gz ${directory}${sample}/${sample}.tmp.fastq.gz

fi

#If present, remove the 3' adapter sequence from the read
#Minimum quality score set to 30 [-q 30]
#Minimum length set to 30 [-m 30]
#Some library preparations may result in sequential ligation of adapters, therefore 2 rounds [-n 2] of adapter removal are performed in case this occured

if [ ! -z "${three_prime_adapter}" ]
then

	cutadapt \
	-a ${three_prime_adapter} \
	-q 30 \
	-m 30 \
	-n 2 \
	-j 12 \
	-o ${directory}${sample}/${sample}.cutadapt.fastq.gz \
	--json=${directory}${sample}/${sample}.cutadapt.3adapter.json \
	${directory}${sample}/${sample}.tmp.fastq.gz

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed

	mv ${directory}${sample}/${sample}.cutadapt.fastq.gz ${directory}${sample}/${sample}.tmp.fastq.gz

fi

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed

	mv ${directory}${sample}/${sample}.tmp.fastq.gz ${directory}${sample}/${sample}.cutadapt.fastq.gz


###Deduplicate Reads###

#Unzip FASTQ file from Cutadapt

	gunzip ${directory}${sample}/${sample}.cutadapt.fastq.gz


##########
#Print the number of reads following nonbarcoded adapter removal in the 'sample.summary.txt' file

	echo "${sample} reads following unbarcoded adapter removal: $(( $(wc -l ${directory}${sample}/${sample}.cutadapt.fastq | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
##########


#Remove duplicate reads, and changes the read ID
#Read ID now contains rank (by duplicate number) and number of duplicates

	awk 'NR%4==2' ${directory}${sample}/${sample}.cutadapt.fastq | \
	sort -T ${directory}${sample} | \
	uniq -c | \
	sort -k1,1nr -T ${directory}${sample} | \
	awk '{print $0,NR}' | \
	awk '{print ">"$3"-"$1"\n"$2}' \
	> ${directory}${sample}/${sample}.cutadapt.deduped.fasta

	rm ${directory}${sample}/${sample}.cutadapt.fastq


##########
#Print the number of reads following ndeduplication in the 'sample.summary.txt' file

	echo "${sample} deduplicated reads: $(( $(wc -l ${directory}${sample}/${sample}.cutadapt.deduped.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
##########


###Remove random barcodes (Cutadapt)###

#If present, remove the 5' barcode sequence from the read
#Minimum quality score set to 30 [-q 30]
#Minimum length set to 30 [-m 30]
#Some library preparations may result in sequential ligation of adapters, therefore 2 rounds [-n 2] of adapter removal are performed in case this occured

if [ ! -z "${five_prime_barcode}" ]
then

	cutadapt \
	-g ^${five_prime_barcode} \
	-m 30 \
	-j 12 \
	-o ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	--json=${directory}${sample}/${sample}.cutadapt.5barcode.json \
	${directory}${sample}/${sample}.cutadapt.deduped.fasta

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed

	mv ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta ${directory}${sample}/${sample}.cutadapt.deduped.fasta

fi

#If present, remove the 3' barcode sequence from the read
#Minimum quality score set to 30 [-q 30]
#Minimum length set to 30 [-m 30]
#Some library preparations may result in sequential ligation of adapters, therefore 2 rounds [-n 2] of adapter removal are performed in case this occured

if [ ! -z "${three_prime_barcode}" ]
then

	cutadapt \
	-a ${three_prime_barcode}$ \
	-m 30 \
	-j 12 \
	-o ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
	--json=${directory}${sample}/${sample}.cutadapt.3barcode.json \
	${directory}${sample}/${sample}.cutadapt.deduped.fasta

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed

	mv ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta ${directory}${sample}/${sample}.cutadapt.deduped.fasta

fi

#Rename the output from Cutadapt so pipeline may proceed, regardless of whether the step was performed

	mv ${directory}${sample}/${sample}.cutadapt.deduped.fasta ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta


##########
#Print the number of reads following barcoded adapter removal in the 'sample.summary.txt' file

	echo "${sample} reads following barcode removal: $(( $(wc -l ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
##########

### Filter our pre-miRNAs and tRNAs (BLAST) ###

if [ $pre_filtered == "yes" ]
then
    echo "Starting pre-miRNA filtering for sample: ${sample}"

    # Align reads to pre-miRNA reference file
    echo "Aligning reads to pre-miRNA reference file..."
    blastn \
        -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
        -query ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
        -out ${directory}${sample}/${sample}.preprocessed.premiRNA.blast \
        -word_size 11 \
        -outfmt 6 \
        -num_threads 1 \
        -strand plus

    echo "BLAST alignment for pre-miRNA completed."

    # Filter pre-miRNA alignments
    echo "Filtering pre-miRNA alignments..."
    awk -v var1=$minimum_pre_miRNA_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.preprocessed.premiRNA.blast | \
        awk -v var1=$minimum_pre_miRNA_bit_score '$12 >= var1' | \
        awk '!x[$1]++' \
        > ${directory}${sample}/${sample}.preprocessed.premiRNA.blast.filtered

    echo "Pre-miRNA filtering completed."

    rm ${directory}${sample}/${sample}.preprocessed.premiRNA.blast

    # Remove reads that aligned to pre-miRNAs with minimum filters
    echo "Removing reads that aligned to pre-miRNAs with minimum filters..."
    join -j1 -v2 \
        <(awk '{print $1}' ${directory}${sample}/${sample}.preprocessed.premiRNA.blast.filtered | sort -k1,1 -T ${directory}${sample} | uniq) \
        <(awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta | sort -k1,1 -T ${directory}${sample}) | \
        awk '{print ">"$1"\n"$2}' \
        > ${directory}${sample}/${sample}.premiRNA.removed.fasta

    ##########
    # Print the number of reads following pre-miRNA filtering in the 'sample.summary.txt' file
    num_reads_premiRNA=$(( $(wc -l ${directory}${sample}/${sample}.premiRNA.removed.fasta | awk '{print $1}') / 2 ))
    echo "${sample} pre-miRNA filtered: ${num_reads_premiRNA}" >> ${directory}${sample}/${sample}.summary.txt
    ##########

    echo "Pre-miRNA removal completed. Number of reads removed: ${num_reads_premiRNA}"

    rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta
    rm ${directory}${sample}/${sample}.preprocessed.premiRNA.blast.filtered

    echo "Aligning reads to tRNA reference file..."
    # Align reads to tRNA reference file
    blastn \
        -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/tRNA_${miRBase_species_abbreviation}.fasta \
        -query ${directory}${sample}/${sample}.premiRNA.removed.fasta \
        -out ${directory}${sample}/${sample}.preprocessed.tRNA.blast \
        -word_size 11 \
        -outfmt 6 \
        -num_threads 1 \
        -strand plus

    echo "BLAST alignment for tRNA completed."

    # Filter tRNA alignments
    echo "Filtering tRNA alignments..."
    awk -v var1=$minimum_tRNA_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.preprocessed.tRNA.blast | \
        awk -v var1=$minimum_tRNA_bit_score '$12 >= var1' | \
        awk '!x[$1]++' \
        > ${directory}${sample}/${sample}.preprocessed.tRNA.blast.filtered

    echo "tRNA filtering completed."

    rm ${directory}${sample}/${sample}.preprocessed.tRNA.blast

    # Remove reads that aligned to tRNAs with minimum filters
    echo "Removing reads that aligned to tRNAs with minimum filters..."
    join -j1 -v2 \
        <(awk '{print $1}' ${directory}${sample}/${sample}.preprocessed.tRNA.blast.filtered | sort -k1,1 -T ${directory}${sample} | uniq) \
        <(awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' ${directory}${sample}/${sample}.premiRNA.removed.fasta | sort -k1,1 -T ${directory}${sample}) | \
        awk '{print ">"$1"\n"$2}' \
        > ${directory}${sample}/${sample}.tRNA.removed.fasta

    ##########
    # Print the number of reads following tRNA filtering in the 'sample.summary.txt' file
    num_reads_tRNA=$(( $(wc -l ${directory}${sample}/${sample}.tRNA.removed.fasta | awk '{print $1}') / 2 ))
    echo "${sample} tRNA filtered: ${num_reads_tRNA}" >> ${directory}${sample}/${sample}.summary.txt
    ##########

    echo "tRNA removal completed. Number of reads removed: ${num_reads_tRNA}"

    rm ${directory}${sample}/${sample}.premiRNA.removed.fasta
    rm ${directory}${sample}/${sample}.preprocessed.tRNA.blast.filtered

    # Rename the output of filtering so that it matches whether or not filtering was performed
    echo "Renaming filtered output to match original file name..."
    mv ${directory}${sample}/${sample}.tRNA.removed.fasta ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta

    echo "Filtering completed for sample: ${sample}"

fi

echo "Starting BLAST alignment to sncRNA reference..."
blastn \
    -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/sncRNA_${miRBase_species_abbreviation}.fasta \
    -query ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
    -out ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
    -word_size 11 \
    -outfmt 6 \
    -num_threads 12

echo "Renaming sncRNA BLAST output..."
mv ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.sncRNA.blast

echo "Running BLAST for miRNA sequences (plus strand)..."
blastn \
    -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
    -query ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
    -out ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
    -word_size 11 \
    -outfmt 6 \
    -num_threads 12 \
    -strand plus

echo "Renaming miRNA plus strand BLAST output..."
mv ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.plus.blast

echo "Running BLAST for miRNA sequences (minus strand)..."
blastn \
    -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
    -query ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
    -out ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
    -word_size 11 \
    -outfmt 6 \
    -num_threads 12 \
    -strand minus

echo "Renaming miRNA minus strand BLAST output..."
mv ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.minus.blast

echo "Running BLAST for hairpin sequences..."
blastn \
    -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
    -query ${reference_directory}/fasta/${miRBase_species_abbreviation}/miRNA_${miRBase_species_abbreviation}.fasta \
    -out ${directory}${sample}/${sample}.hairpin.miRNA.blast \
    -word_size 11 \
    -outfmt 6 \
    -num_threads 12

echo "Running BLAST for hairpin sequences (plus strand)..."
blastn \
    -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
    -query ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
    -out ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
    -word_size 11 \
    -outfmt 6 \
    -num_threads 12 \
    -strand plus

echo "Renaming hairpin plus strand BLAST output..."
mv ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.plus.blast

echo "Running BLAST for hairpin sequences (minus strand)..."
blastn \
    -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
    -query ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
    -out ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
    -word_size 11 \
    -outfmt 6 \
    -num_threads 12 \
    -strand minus

echo "Renaming hairpin minus strand BLAST output..."
mv ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.minus.blast

echo "Filtering BLAST results..."
awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.plus.blast | \
    awk -v var1=$minimum_length '$4 >= var1' | \
    awk -v var1=$maximum_mismatch '$5 <= var1' | \
    awk -v var1=$maximum_gap '$6 <= var1' | \
    awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
    awk '!x[$1]++' \
    > ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.plus.blast.filtered

echo "Filtering BLAST results for miRNA minus strand..."
awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.minus.blast | \
    awk -v var1=$minimum_length '$4 >= var1' | \
    awk -v var1=$maximum_mismatch '$5 <= var1' | \
    awk -v var1=$maximum_gap '$6 <= var1' | \
    awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
    awk '!x[$1]++' \
    > ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.minus.blast.filtered

echo "Filtering BLAST results for hairpin plus strand..."
awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.plus.blast | \
    awk -v var1=$minimum_length '$4 >= var1' | \
    awk -v var1=$maximum_mismatch '$5 <= var1' | \
    awk -v var1=$maximum_gap '$6 <= var1' | \
    awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
    awk '!x[$1]++' \
    > ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.plus.blast.filtered

echo "Filtering BLAST results for hairpin minus strand..."
awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.minus.blast | \
    awk -v var1=$minimum_length '$4 >= var1' | \
    awk -v var1=$maximum_mismatch '$5 <= var1' | \
    awk -v var1=$maximum_gap '$6 <= var1' | \
    awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
    awk '!x[$1]++' \
    > ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.minus.blast.filtered

echo "Merging filtered results..."
join \
    <(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.plus.blast.filtered) \
    <(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.plus.blast.filtered) \
    > ${directory}${sample}/${sample}.blast.miRNA_iso.plus.merged

join \
    <(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.miRNA.minus.blast.filtered) \
    <(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.hairpin.minus.blast.filtered) \
    > ${directory}${sample}/${sample}.blast.miRNA_iso.minus.merged

echo "BLAST pipeline completed successfully!"

####################################################################################
# The objective is to get the isomiRic miRNA alignments

# Join miRNA alignments with hairpin miRNA alignments
echo "Joining miRNA alignments with hairpin alignments (plus strand)..."
join -1 13 -2 2 \
    <(sort -k13 ${directory}${sample}/${sample}.blast.miRNA_iso.plus.merged) \
    <(sort -k2 ${directory}${sample}/${sample}.hairpin.miRNA.blast) \
    | awk '{if($3==$24)print $0}' | sort -k2 > temp1.plus
echo "Joining miRNA alignments with hairpin alignments (minus strand)..."
join -1 13 -2 2 \
    <(sort -k13 ${directory}${sample}/${sample}.blast.miRNA_iso.minus.merged) \
    <(sort -k2 ${directory}${sample}/${sample}.hairpin.miRNA.blast) \
    | awk '{if($3==$24)print $0}' | sort -k2 > temp1.minus

# Check if the alignment length is different
echo "Checking alignment length and direction for plus strand..."
awk 'function abs(x){return ((x < 0.0) ? -x : x)} \
		{ \
			if($26 >= 17){ \
				if($32 > $31){ \
					if($32 > $21 && $31 > $20) \
						print $0,"\t3p:-"abs($32-$21)";5p:+"abs($31-$20); \
					else if($32 == $21 && $31 > $20) \
						print $0,"\t5p:+"abs($31-$20); \
					else if($32 == $21 && $31 < $20) \
						print $0,"\t5p:-"abs($31-$20); \
					else if($32 < $21 && $31 > $20) \
						print $0 ,"\t3p:+"abs($32-$21)";5p:+"abs($31-$20); \
					else if($32 < $21 && $31 == $20) \
						print $0 ,"\t3p:+"abs($32-$21); \
					else if($32 > $21 && $31 < $20) \
						print $0,"\t3p:-"abs($32-$21)";5p:-"abs($31-$20); \
					else if($32 > $21 && $31 == $20) \
						print $0,"\t3p:-"abs($32-$21); \
					else if($32 < $21 && $31 < $20) \
						print $0 ,"\t3p:+"abs($32-$21)";5p:-"abs($31-$20); \
					else \
						print $0; \
				} \
			} else \
				print $0"\tND\tmiRNA alignment length"; \
		}' temp1.plus > plus
echo "Checking alignment length and direction for minus strand..."
awk 'function abs(x){return ((x < 0.0) ? -x : x)} \
		{ \
			if($26 >= 17){ \
				if($32 > $31){ \
					if($32 > $20 && $31 > $21) \
						print $0,"\t3p:-"abs($32-$20)";5p:+"abs($31-$21); \
					else if($32 == $20 && $31 > $21) \
						print $0,"\t5p:+"abs($31-$21); \
					else if($32 == $20 && $31 < $21) \
						print $0,"\t5p:-"abs($31-$21); \
					else if($32 < $20 && $31 > $21) \
						print $0 ,"\t3p:+"abs($32-$20)";5p:+"abs($31-$21); \
					else if($32 < $20 && $31 == $21) \
						print $0 ,"\t3p:+"abs($32-$20); \
					else if($32 > $20 && $31 < $21) \
						print $0,"\t3p:-"abs($32-$20)";5p:-"abs($31-$21); \
					else if($32 > $20 && $31 == $21) \
						print $0,"\t3p:-"abs($32-$20); \
					else if($32 < $20 && $31 < $21) \
						print $0 ,"\t3p:+"abs($32-$20)";5p:-"abs($31-$21); \
					else \
						print $0; \
				} \
			} else \
				print $0"\tND\tmiRNA alignment length"; \
		}' temp1.minus > minus

# Pull out the canonical miRNA alignments from both files, and make the alignment file for further processing
echo "Extracting canonical miRNA alignments for plus strand..."
awk '{if(NF==34)print $2"\t"$3"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' plus > plus.miRNA
echo "Extracting canonical miRNA alignments for minus strand..."
awk '{if(NF==34)print $2"\t"$3"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' minus > minus.miRNA

# Pull out the isomiRic miRNA alignments from both files
echo "Extracting isomiRic miRNA alignments for plus strand..."
awk '{if(NF==35)print $2"\t"$3"|"$NF"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' plus > plus.miRNA_isomiR
echo "Extracting isomiRic miRNA alignments for minus strand..."
awk '{if(NF==35)print $2"\t"$3"|"$NF"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23}' minus > minus.miRNA_isomiR

# Put the resulting alignment in the total annotation file
echo "Combining all miRNA and isomiR alignments into one file..."
cat plus.miRNA minus.miRNA plus.miRNA_isomiR minus.miRNA_isomiR > test.miRNA.alignments.blast

# Make tabular version of barcoded FASTA file (from before BLAST)
echo "Converting barcoded FASTA to tabular format..."
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' \
    ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
    > ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab

# Rename the file to join the pipeline
echo "Renaming the final alignment file..."
mv test.miRNA.alignments.blast ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered

echo "Processing completed!"
####################################################################################

# Join tabular input file with filtered BLAST results
echo "Joining BLAST filtered results with tabular data for sample: ${sample}"

join \
	<(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast.filtered -T ${directory}${sample}) \
	<(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab -T ${directory}${sample}) \
	> ${directory}${sample}/${sample}.blast.merged

echo "Joining completed for sample: ${sample}"

##########
# Print the number of reads containing a sncRNA (satisfying BLAST filters) in the 'sample.summary.txt' file
echo "Counting filtered BLAST alignments for sample: ${sample}"

echo "${sample} filtered BLAST alignments: $(( $(wc -l ${directory}${sample}/${sample}.blast.merged | awk '{print $1}') / 1 ))" >> ${directory}${sample}/${sample}.summary.txt

echo "Filtered BLAST alignment count written to summary file"
##########

echo "Proceeding to identify potential chimeras for sample: ${sample}"

### Identify potential chimeras ###
# Create FASTA file with sequence of sncRNA removed and name appended to the read ID
# Only output reads containing minimum length of sequence following sncRNA removal

awk '{stop=$8;readlength=length($13);print $0,"\t", readlength-(stop)}' ${directory}${sample}/${sample}.blast.merged | \
	awk -v var1=$minimum_length_after_sncRNA '$14 >= var1' | \
	awk '{print $0,"\t",substr($13, $8+1, length($13)-($8))}' | \
	awk '{print">"$1"."$2"\n"$15}' \
	> ${directory}${sample}/${sample}.target.fasta

echo "FASTA file for potential chimeras created"

##########
# Print the number of reads containing a sncRNA (satisfying BLAST filters) and minimum sequence length
echo "Counting potential chimeras for sample: ${sample}"

echo "${sample} potential chimeras: $(( $(wc -l ${directory}${sample}/${sample}.target.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt

echo "Potential chimeras count written to summary file"
##########

echo "Proceeding to alignment step for sample: ${sample}"

### Align reads to the reference genome (HISAT2) ###
echo "Running HISAT2 alignment for sample: ${sample}"

hisat2 \
	-x ${reference_directory}/fasta/${genome_species_abbreviation}/${genome_species_abbreviation} \
	-f ${directory}${sample}/${sample}.target.fasta \
	-S ${directory}${sample}/${sample}.aligned.sam \
	--summary-file ${directory}${sample}/${sample}.hisat2summary.txt

echo "HISAT2 alignment completed for sample: ${sample}"

##########
# Print the number of reads that aligned to the genome a single time
echo "Counting unique alignments for sample: ${sample}"

echo "${sample} unique alignments: $(( $(awk '$5 == 60' ${directory}${sample}/${sample}.aligned.sam | wc -l | awk '{print $1}') / 1 ))" >> ${directory}${sample}/${sample}.summary.txt

echo "Unique alignment count written to summary file"
##########

echo "Filtering uniquely mapped reads for sample: ${sample}"

awk '/^@/ || $5 == 60' ${directory}${sample}/${sample}.aligned.sam \
	> ${directory}${sample}/${sample}.aligned.unique.sam

echo "Unique reads filtered"

# Convert to BAM file
echo "Converting SAM to BAM format for sample: ${sample}"

samtools view \
	-S \
	-h \
	-b ${directory}${sample}/${sample}.aligned.unique.sam \
	> ${directory}${sample}/${sample}.aligned.unique.bam

echo "BAM file created for sample: ${sample}"

echo "Finish: $(date)" >> ${directory}${sample}/${sample}.summary.txt
echo "${sample} completed successfully"

done

echo "All samples processed. Deactivating conda environment."
conda deactivate
