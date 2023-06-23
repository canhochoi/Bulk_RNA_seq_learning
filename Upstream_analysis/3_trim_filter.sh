#!/usr/bin/bash

p_sra='/home/lqluan/bulk_RNAseq/week7/sacCer2/sra'
p_raw='/home/lqluan/bulk_RNAseq/week7/sacCer2/raw'
p_trim='/home/lqluan/bulk_RNAseq/week7/sacCer2/trim'
trimmomatic_jar='/home/lqluan/bulk_RNAseq/week7/tools/Trimmomatic-0.39/trimmomatic-0.39.jar'

for file in $(ls *_R1*.gz);do
    # Extract the base name without the extension
    base_name="${file%_*.*.*}"

    # Construct the input and output file names
    read_1=${base_name}_R1.fastq.gz
    read_2=${base_name}_R2.fastq.gz
    output_paired_1="${base_name}_R1_paired.fastq.gz"
    output_unpaired_1="${base_name}_R1_unpaired.fastq.gz"
    output_paired_2="${base_name}_R2_paired.fastq.gz"
    output_unpaired_2="${base_name}_R2_unpaired.fastq.gz"
    
    # Run Trimmomatic with the specified parameters
    java -jar $trimmomatic_jar PE \
        -threads 4 \
        -trimlog "$p_trim/trim_test.log" \
        "$p_raw/$read_1" \
        "$p_raw/$read_2" \
        "$p_trim/$output_paired_1" \
        "$p_trim/$output_unpaired_1" \
        "$p_trim/$output_paired_2" \
        "$p_trim/$output_unpaired_2" \
        HEADCROP:10 \
        CROP:65 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:30 \
        MINLEN:36
    
    # Remove unpaired file
    rm $p_trim/$output_unpaired_1 $p_trim/$output_unpaired_2

    # Check fastqc again
    fastqc $p_trim/$output_paired_1 -o $p_trim/qc_check
    fastqc $p_trim/$output_paired_2 -o $p_trim/qc_check
done
