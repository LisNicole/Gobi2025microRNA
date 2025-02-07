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



###Adapter Extraction###

#Getting sample names from adapter file
samples=$(awk '!/^#/' ${adapter_file} | awk '{print $1}')


for sample in $samples
do

  #Extract adapter and barcode sequences from adapter file
five_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $2}')
three_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $3}')
five_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $4}')
three_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $5}')

  #Remove 'sample.summary.txt' file if it exists, to ensure a proper end file

    [ -e ${directory}${sample}/${sample}.summary.txt ] && rm ${directory}${sample}/${sample}.summary.txt

  #Print the name of the current sample being processed in the 'sample.summary.txt' file

    echo "$sample" >> ${directory}${sample}/${sample}.summary.txt

  #Print the structure of the reads for the sample being processed in the 'sample.summary.txt' file

    echo "5'-${five_prime_adapter}${five_prime_barcode}...miRNA...${three_prime_barcode}${three_prime_adapter}-3'" >> ${directory}${sample}/${sample}.summary.txt


  #Print the date and time that processing begins on the sample being processed in the 'sample.summary.txt' file

    echo "Start: $(date)" >> ${directory}${sample}/${sample}.summary.txt


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
#IF the user set the paired-end parameter [-p] to anything other than yes skip this step

#If there are paired end reads
if [ $paired_end == "yes" ]
then

#Make directory for FLASh output

	mkdir ${directory}${sample}/${sample}_FLASh

#Run FLASh

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


###Cutadapt to Remove Adapters without Random Nucleotides, Quality Trimming and Minimum Threshold###

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
#Print the number of reads following deduplication in the 'sample.summary.txt' file

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

done ####temporary####