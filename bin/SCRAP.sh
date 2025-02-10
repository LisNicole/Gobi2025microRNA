#!/bin/sh

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

# Add trailing slash if missing
directory=$(echo "${directory}" | sed 's![^/]$!&/!')
reference_directory=$(echo "${reference_directory}" | sed 's![^/]$!&/!')

##############################
# PARAMETERS
##############################

# Parameters for filtering pre-miRNA BLAST results
minimum_pre_miRNA_evalue=.05
minimum_pre_miRNA_bit_score=50

# Parameters for filtering tRNA BLAST results
minimum_tRNA_evalue=.05
minimum_tRNA_bit_score=50

# Parameters for filtering rRNA BLAST results
minimum_rRNA_evalue=.05
minimum_rRNA_bit_score=50

# Parameters for filtering primer BLAST results
minimum_primer_evalue=.05
minimum_primer_bit_score=50

# ---------------------------
# PARAMETERS for filtering BLAST results (sncRNA/isomiR detection)
# Modified for isomiR detection:
minimum_evalue=.05
minimum_length=14         # Changed from 20 to 14
maximum_mismatch=2        # Changed from 1 to 2
maximum_gap=1
maximum_gap_if_maximum_mismatch=0

# Minimum length of sequence after the sncRNA within a read (after adapter and barcode removal)
minimum_length_after_sncRNA=15

##############################################################################
# The remainder of the pipeline (QC, adapter removal, deduplication, etc.) remains as in the original script
##############################################################################

# Activate conda environment
location=$(conda info | awk '/base environment/' | awk '{print $4}')
source ${location}/etc/profile.d/conda.sh
conda activate SCRAP

samples=$(awk '!/^#/' ${adapter_file} | awk '{print $1}')

### Quality Control (FastQC) ###
mkdir ${directory}FastQC_Reports
for sample in $samples
do
    for file in $(ls ${directory}${sample}/ | awk '/fastq/')
    do
        fastqc -o ${directory}FastQC_Reports ${directory}${sample}/${file}
    done
done

### Quality Control (MultiQC) ###
mkdir ${directory}MultiQC_Report
multiqc ${directory}FastQC_Reports -o ${directory}MultiQC_Report

for sample in $samples
do
    # Extract adapter and barcode sequences from adapter file
    five_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $2}')
    three_prime_adapter=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $3}')
    five_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $4}')
    three_prime_barcode=$(awk -v var1=$sample '$1==var1' ${adapter_file} | awk -F '\t' '{print $5}')

    # Remove sample summary file if it exists
    [ -e ${directory}${sample}/${sample}.summary.txt ] && rm ${directory}${sample}/${sample}.summary.txt

    echo "$sample" >> ${directory}${sample}/${sample}.summary.txt
    echo "5'-${five_prime_adapter}${five_prime_barcode}...sncRNA-targetRNA...${three_prime_barcode}${three_prime_adapter}-3'" >> ${directory}${sample}/${sample}.summary.txt
    echo "miRBase Species Abbreviation: ${miRBase_species_abbreviation}" >> ${directory}${sample}/${sample}.summary.txt
    echo "Genome Species Abbreviation: ${genome_species_abbreviation}" >> ${directory}${sample}/${sample}.summary.txt
    echo "Start: $(date)" >> ${directory}${sample}/${sample}.summary.txt

    # Count raw reads for each file
    for file in $(ls ${directory}${sample}/ | awk '/fastq/')
    do
        gunzip -c ${directory}${sample}/${file} > ${directory}${sample}/${file}.tmp
        echo "${file} raw reads: $(( $(wc -l ${directory}${sample}/${file}.tmp | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
        rm ${directory}${sample}/${file}.tmp
    done

    ### Combine Paired-End Reads (FLASh) ###
    if [ $paired_end == "yes" ]
    then
        mkdir ${directory}${sample}/${sample}_FLASh
        flash --allow-outies \
              --output-directory=${directory}${sample}/${sample}_FLASh/ \
              --output-prefix=${sample} \
              --max-overlap=150 \
              --min-overlap=6 \
              --compress \
              ${directory}${sample}/${sample}_R1.fastq.gz \
              ${directory}${sample}/${sample}_R2.fastq.gz \
              2>&1 | tee ${directory}${sample}/${sample}_FLASh/FLASh_${sample}.log

        mv ${directory}${sample}/${sample}_FLASh/${sample}.extendedFrags.fastq.gz ${directory}${sample}/${sample}.fastq.gz
        [ -e ${directory}${sample}/${sample}_FLASh ] && rm -r ${directory}${sample}/${sample}_FLASh

        gunzip -k ${directory}${sample}/${sample}.fastq.gz
        echo "${sample} combined paired-end reads: $(( $(wc -l ${directory}${sample}/${sample}.fastq | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt
        rm ${directory}${sample}/${sample}.fastq
    fi

    cp ${directory}${sample}/${sample}.fastq.gz ${directory}${sample}/${sample}.tmp.fastq.gz

    ### Remove Adapters, Quality Trimming, and Minimum Length Filtering (Cutadapt) ###
    if [ ! -z "${five_prime_adapter}" ]
    then
        cutadapt -g ${five_prime_adapter} -q 30 -m 30 -n 2 \
                 -o ${directory}${sample}/${sample}.cutadapt.fastq.gz \
                 --json=${directory}${sample}/${sample}.cutadapt.5adapter.json \
                 ${directory}${sample}/${sample}.tmp.fastq.gz
        mv ${directory}${sample}/${sample}.cutadapt.fastq.gz ${directory}${sample}/${sample}.tmp.fastq.gz
    fi

    if [ ! -z "${three_prime_adapter}" ]
    then
        cutadapt -a ${three_prime_adapter} -q 30 -m 30 -n 2 \
                 -o ${directory}${sample}/${sample}.cutadapt.fastq.gz \
                 --json=${directory}${sample}/${sample}.cutadapt.3adapter.json \
                 ${directory}${sample}/${sample}.tmp.fastq.gz
        mv ${directory}${sample}/${sample}.cutadapt.fastq.gz ${directory}${sample}/${sample}.tmp.fastq.gz
    fi

    mv ${directory}${sample}/${sample}.tmp.fastq.gz ${directory}${sample}/${sample}.cutadapt.fastq.gz

    ### Deduplicate Reads ###
    gunzip ${directory}${sample}/${sample}.cutadapt.fastq.gz
    echo "${sample} reads following unbarcoded adapter removal: $(( $(wc -l ${directory}${sample}/${sample}.cutadapt.fastq | awk '{print $1}') / 4 ))" >> ${directory}${sample}/${sample}.summary.txt

    awk 'NR%4==2' ${directory}${sample}/${sample}.cutadapt.fastq | \
      sort -T ${directory}${sample} | uniq -c | sort -k1,1nr -T ${directory}${sample} | \
      awk '{print $0,NR}' | \
      awk '{print ">"$3"-"$1"\n"$2}' \
      > ${directory}${sample}/${sample}.cutadapt.deduped.fasta
    rm ${directory}${sample}/${sample}.cutadapt.fastq

    echo "${sample} deduplicated reads: $(( $(wc -l ${directory}${sample}/${sample}.cutadapt.deduped.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt

    ### Remove random barcodes (Cutadapt) ###
    if [ ! -z "${five_prime_barcode}" ]
    then
        cutadapt -g ^${five_prime_barcode} -m 30 \
                 -o ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
                 --json=${directory}${sample}/${sample}.cutadapt.5barcode.json \
                 ${directory}${sample}/${sample}.cutadapt.deduped.fasta
        mv ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta ${directory}${sample}/${sample}.cutadapt.deduped.fasta
    fi

    if [ ! -z "${three_prime_barcode}" ]
    then
        cutadapt -a ${three_prime_barcode}$ -m 30 \
                 -o ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
                 --json=${directory}${sample}/${sample}.cutadapt.3barcode.json \
                 ${directory}${sample}/${sample}.cutadapt.deduped.fasta
        mv ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta ${directory}${sample}/${sample}.cutadapt.deduped.fasta
    fi

    mv ${directory}${sample}/${sample}.cutadapt.deduped.fasta ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta

    echo "${sample} reads following barcode removal: $(( $(wc -l ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt

    ################### SEARCH FOR sncRNAs ################################

    ### Identify sncRNA (BLAST) ###
    blastn -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/sncRNA_${miRBase_species_abbreviation}.fasta \
           -query ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
           -out ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast \
           -word_size 11 \
           -outfmt 6 \
           -num_threads 1 \
           -strand plus

    # Filter BLAST results using the modified parameters
    awk -v var1=$minimum_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast | \
      awk -v var1=$minimum_length '$4 >= var1' | \
      awk -v var1=$maximum_mismatch '$5 <= var1' | \
      awk -v var1=$maximum_gap '$6 <= var1' | \
      awk -v var1=$maximum_mismatch -v var2=$maximum_gap_if_maximum_mismatch '!($5 == var1 && $6 > var2)' | \
      awk '!x[$1]++' \
      > ${directory}${sample}/${sample}.sncrnaidentified.blast

    rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.blast

    # Create a tab-delimited version of the FASTA file (for joining with BLAST results)
    awk 'BEGIN{RS=">";OFS="\t"} NR>1{print $1,$2}' \
         ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta \
         > ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab
    rm ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.fasta

    # Join the BLAST output with the FASTA records
    join <(sort -k1,1 ${directory}${sample}/${sample}.sncrnaidentified.blast -T ${directory}${sample}) \
         <(sort -k1,1 ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.tab -T ${directory}${sample}) \
         > ${directory}${sample}/${sample}.blast.merged

    echo "${sample} detected sncRNAs: $(( $(wc -l ${directory}${sample}/${sample}.blast.merged | awk '{print $1}') ))" >> ${directory}${sample}/${sample}.summary.txt

    ############################################################################
    # NEW: ISOMIR EXTRACTION STEP
    #
    # Two new outputs are produced:
    # (A) A tab-delimited file containing:
    #     1. Read ID
    #     2. sncRNA (miRNA) ID
    #     3. Query start (read coordinate of alignment start)
    #     4. Query end (read coordinate of alignment end)
    #     5. Subject start (miRNA reference coordinate of alignment start)
    #     6. Subject end (miRNA reference coordinate of alignment end)
    #     7. Number of mismatches
    #     8. Number of gap opens
    #     9. Extracted isomiR sequence (from read, positions [qstart, qend])
    #    10. Target sequence (read sequence after the sncRNA)
    #    11. Full read length
    #
    # (B) A FASTA file with a header that includes alignment info and a sequence
    #     containing the extracted isomiR.
    ############################################################################

    awk 'BEGIN{OFS="\t"} {
         # Fields in the blast merged file (assuming they are in the same order as SCRAP’s original):
         # $1 = Read ID, $2 = sncRNA ID, $4 = alignment length, $5 = mismatches, $6 = gap count,
         # $7 = query start, $8 = query end, $9 = subject start, $10 = subject end,
         # $13 = full read sequence.
         isomir_seq = substr($13, $7, $8 - $7 + 1);
         target_seq = substr($13, $8+1);
         read_len = length($13);
         print $1, $2, $7, $8, $9, $10, $5, $6, isomir_seq, target_seq, read_len
    }' ${directory}${sample}/${sample}.blast.merged \
      > ${directory}${sample}/${sample}.isomir_details.txt

    awk '{
         isomir = substr($13, $7, $8 - $7 + 1);
         print ">"$1"."$2"|qstart:"$7"|qend:"$8"|mm:"$5"|gaps:"$6"\n" isomir
    }' ${directory}${sample}/${sample}.blast.merged \
      > ${directory}${sample}/${sample}.isomir.fasta

    ### Identify potential chimeras ###
    # Create FASTA file with the sequence from the beginning of the sncRNA (in the read) to the end of the read.
    awk '{start=$7; readlength=length($13); print $0,"\t", readlength-(start)}' ${directory}${sample}/${sample}.blast.merged | \
      awk -v var1=$minimum_length_after_sncRNA '$14 >= var1' | \
      awk '{print $0,"\t", substr($13, $7, length($13)-($7)+1)}' | \
      awk '{print ">"$1"."$2"\n"$15}' \
      > ${directory}${sample}/${sample}.wholechimera.fasta

    # Create FASTA file with the sequence after the sncRNA (target only)
    awk '{stop=$8; readlength=length($13); print $0,"\t", readlength-(stop)}' ${directory}${sample}/${sample}.blast.merged | \
      awk -v var1=$minimum_length_after_sncRNA '$14 >= var1' | \
      awk '{print $0,"\t", substr($13, $8+1, length($13)-($8))}' | \
      awk '{print ">"$1"."$2"\n"$15}' \
      > ${directory}${sample}/${sample}.target.fasta

    rm ${directory}${sample}/${sample}.blast.merged

    echo "${sample} potential chimeras: $(( $(wc -l ${directory}${sample}/${sample}.wholechimera.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt

    ### Filter out pre-miRNAs and tRNAs (BLAST) ###
    if [ $pre_filtered == "yes" ]
    then
        # Pre-miRNA filtering
        blastn -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/hairpin_${miRBase_species_abbreviation}.fasta \
               -query ${directory}${sample}/${sample}.wholechimera.fasta \
               -out ${directory}${sample}/${sample}.preprocessed.premiRNA.blast \
               -word_size 11 \
               -outfmt 6 \
               -num_threads 1 \
               -strand plus

        awk -v var1=$minimum_pre_miRNA_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.preprocessed.premiRNA.blast | \
          awk -v var1=$minimum_pre_miRNA_bit_score '$12 >= var1' | \
          awk '$7 < 14' | \
          awk '!x[$1]++' \
          > ${directory}${sample}/${sample}.preprocessed.premiRNA.blast.filtered

        rm ${directory}${sample}/${sample}.preprocessed.premiRNA.blast

        join -j1 -v2 <(awk '{print $1}' ${directory}${sample}/${sample}.preprocessed.premiRNA.blast.filtered | sort -k1,1 -T ${directory}${sample} | uniq) \
             <(awk 'BEGIN{RS=">";OFS="\t"} NR>1{print $1,$2}' ${directory}${sample}/${sample}.wholechimera.fasta | sort -k1,1 -T ${directory}${sample}) | \
             awk '{print ">"$1"\n"$2}' \
             > ${directory}${sample}/${sample}.premiRNA.removed.fasta

        echo "${sample} pre-miRNA filtered: $(( $(wc -l ${directory}${sample}/${sample}.premiRNA.removed.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
        rm ${directory}${sample}/${sample}.preprocessed.premiRNA.blast.filtered

        # tRNA filtering
        blastn -db ${reference_directory}/fasta/${miRBase_species_abbreviation}/tRNA_${miRBase_species_abbreviation}.fasta \
               -query ${directory}${sample}/${sample}.premiRNA.removed.fasta \
               -out ${directory}${sample}/${sample}.preprocessed.tRNA.blast \
               -word_size 11 \
               -outfmt 6 \
               -num_threads 1 \
               -strand plus

        awk -v var1=$minimum_tRNA_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.preprocessed.tRNA.blast | \
          awk -v var1=$minimum_tRNA_bit_score '$12 >= var1' | \
          awk '$7 < 14' | \
          awk '!x[$1]++' \
          > ${directory}${sample}/${sample}.preprocessed.tRNA.blast.filtered

        rm ${directory}${sample}/${sample}.preprocessed.tRNA.blast

        join -j1 -v2 <(awk '{print $1}' ${directory}${sample}/${sample}.preprocessed.tRNA.blast.filtered | sort -k1,1 -T ${directory}${sample} | uniq) \
             <(awk 'BEGIN{RS=">";OFS="\t"} NR>1{print $1,$2}' ${directory}${sample}/${sample}.premiRNA.removed.fasta | sort -k1,1 -T ${directory}${sample}) | \
             awk '{print ">"$1"\n"$2}' \
             > ${directory}${sample}/${sample}.tRNA.removed.fasta

        echo "${sample} tRNA filtered: $(( $(wc -l ${directory}${sample}/${sample}.tRNA.removed.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
        rm ${directory}${sample}/${sample}.premiRNA.removed.fasta
        rm ${directory}${sample}/${sample}.preprocessed.tRNA.blast.filtered

        # rRNA filtering
        blastn -db ${reference_directory}/test_filter/mouse_rRNA.fasta \
               -query ${directory}${sample}/${sample}.tRNA.removed.fasta \
               -out ${directory}${sample}/${sample}.preprocessed.rRNA.blast \
               -word_size 11 \
               -outfmt 6 \
               -num_threads 1 \
               -strand plus

        awk -v var1=$minimum_rRNA_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.preprocessed.rRNA.blast | \
          awk -v var1=$minimum_rRNA_bit_score '$12 >= var1' | \
          awk '$7 < 14' | \
          awk '!x[$1]++' \
          > ${directory}${sample}/${sample}.preprocessed.rRNA.blast.filtered

        rm ${directory}${sample}/${sample}.preprocessed.rRNA.blast

        join -j1 -v2 <(awk '{print $1}' ${directory}${sample}/${sample}.preprocessed.rRNA.blast.filtered | sort -k1,1 -T ${directory}${sample} | uniq) \
             <(awk 'BEGIN{RS=">";OFS="\t"} NR>1{print $1,$2}' ${directory}${sample}/${sample}.tRNA.removed.fasta | sort -k1,1 -T ${directory}${sample}) | \
             awk '{print ">"$1"\n"$2}' \
             > ${directory}${sample}/${sample}.rRNA.removed.fasta

        echo "${sample} rRNA filtered: $(( $(wc -l ${directory}${sample}/${sample}.rRNA.removed.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
        rm ${directory}${sample}/${sample}.tRNA.removed.fasta
        rm ${directory}${sample}/${sample}.preprocessed.rRNA.blast.filtered

        # Primer filtering
        blastn -db ${reference_directory}/test_filter/CIMERA-seq_primer.fasta \
               -query ${directory}${sample}/${sample}.rRNA.removed.fasta \
               -out ${directory}${sample}/${sample}.preprocessed.primer.blast \
               -word_size 11 \
               -outfmt 6 \
               -num_threads 1 \
               -strand plus

        awk -v var1=$minimum_primer_evalue '($11 + 0) <= var1' ${directory}${sample}/${sample}.preprocessed.primer.blast | \
          awk -v var1=$minimum_primer_bit_score '$12 >= var1' | \
          awk '$7 < 14' | \
          awk '!x[$1]++' \
          > ${directory}${sample}/${sample}.preprocessed.primer.blast.filtered

        rm ${directory}${sample}/${sample}.preprocessed.primer.blast

        join -j1 -v2 <(awk '{print $1}' ${directory}${sample}/${sample}.preprocessed.primer.blast.filtered | sort -k1,1 -T ${directory}${sample} | uniq) \
             <(awk 'BEGIN{RS=">";OFS="\t"} NR>1{print $1,$2}' ${directory}${sample}/${sample}.rRNA.removed.fasta | sort -k1,1 -T ${directory}${sample}) | \
             awk '{print ">"$1"\n"$2}' \
             > ${directory}${sample}/${sample}.primer.removed.fasta

        echo "${sample} primer filtered: $(( $(wc -l ${directory}${sample}/${sample}.primer.removed.fasta | awk '{print $1}') / 2 ))" >> ${directory}${sample}/${sample}.summary.txt
        rm ${directory}${sample}/${sample}.rRNA.removed.fasta
        rm ${directory}${sample}/${sample}.preprocessed.primer.blast.filtered

        mv ${directory}${sample}/${sample}.tRNA.removed.fasta ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.filtered.fasta
        mv ${directory}${sample}/${sample}.primer.removed.fasta ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.filtered.fasta
    fi

    # Remove any chimeras that were filtered out from your target file
    awk 'BEGIN{RS=">";OFS="\t"} NR>1{print $1,$2}' ${directory}${sample}/${sample}.target.fasta | sort -k1,1 -T ${directory}${sample} \
         > ${directory}${sample}/${sample}.prejoin.target.fasta
    awk 'BEGIN{RS=">";OFS="\t"} NR>1{print $1,$2}' ${directory}${sample}/${sample}.cutadapt.deduped.barcoded.filtered.fasta | sort -k1,1 -T ${directory}${sample} \
         > ${directory}${sample}/${sample}.prejoin.chimeras.fasta

    join -j1 -o 1.1,1.2 ${directory}${sample}/${sample}.prejoin.target.fasta ${directory}${sample}/${sample}.prejoin.chimeras.fasta | \
      awk '{print ">"$1"\n"$2}' \
      > ${directory}${sample}/${sample}.filtered.target.fasta

    ### Align reads to the reference genome (HISAT2) ###
    hisat2 -x ${reference_directory}/fasta/${genome_species_abbreviation}/${genome_species_abbreviation} \
           -f ${directory}${sample}/${sample}.filtered.target.fasta \
           -S ${directory}${sample}/${sample}.aligned.sam \
           --summary-file ${directory}${sample}/${sample}.hisat2summary.txt

    echo "${sample} unique alignments: $(( $(awk '$5 == 60' ${directory}${sample}/${sample}.aligned.sam | wc -l | awk '{print $1}') ))" >> ${directory}${sample}/${sample}.summary.txt

    awk '/^@/ || $5 == 60' ${directory}${sample}/${sample}.aligned.sam \
         > ${directory}${sample}/${sample}.aligned.unique.sam

    samtools view -S -h -b ${directory}${sample}/${sample}.aligned.unique.sam \
         > ${directory}${sample}/${sample}.aligned.unique.bam

    echo "Finish: $(date)" >> ${directory}${sample}/${sample}.summary.txt

done

conda deactivate
