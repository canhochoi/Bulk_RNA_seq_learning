#!/usr/bin/bash

p_sra='/home/lqluan/bulk_RNAseq/week7/sacCer2/sra'
p_raw='/home/lqluan/bulk_RNAseq/week7/sacCer2/raw'
for sra_acc in $(cat $p_sra/SraAccList.csv); do

 # Get the SRA
 prefetch -v $sra_acc --output-directory $p_sra
    echo "=== $sra_acc sra file downloaded ==="
 # Download fastq file
 fastq-dump \
    --outdir $p_raw/ \
    --split-files $p_sra/${sra_acc}/${sra_acc}.sra \
    && gzip -f $p_raw/*.fastq
    echo "=== Split $sra_acc file completed ==="

done

cd $p_raw
mv SRR23867633_1.fastq.gz WT_C_1_R1.fastq.gz
mv SRR23867633_2.fastq.gz WT_C_1_R2.fastq.gz
mv SRR23867634_1.fastq.gz WT_E_1_R1.fastq.gz
mv SRR23867634_2.fastq.gz WT_E_1_R2.fastq.gz
mv SRR23867635_1.fastq.gz WT_C_2_R1.fastq.gz
mv SRR23867635_2.fastq.gz WT_C_2_R2.fastq.gz
mv SRR23867632_1.fastq.gz WT_E_2_R1.fastq.gz
mv SRR23867632_2.fastq.gz WT_E_2_R2.fastq.gz
mv SRR23867636_1.fastq.gz WT_C_3_R1.fastq.gz
mv SRR23867636_2.fastq.gz WT_C_3_R2.fastq.gz
mv SRR23867637_1.fastq.gz WT_E_3_R1.fastq.gz
mv SRR23867637_2.fastq.gz WT_E_3_R2.fastq.gz


