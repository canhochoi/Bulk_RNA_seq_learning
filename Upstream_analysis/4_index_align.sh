#!/usr/bin/bash
p_ref='/home/lqluan/bulk_RNAseq/week7/sacCer2/ref'
p_align='/home/lqluan/bulk_RNAseq/week7/sacCer2/align'
p_trim='/home/lqluan/bulk_RNAseq/week7/sacCer2/trim'

# Using STAR to index the genome
STAR \
--runThreadN 12 \
--runMode genomeGenerate \
--genomeDir $p_ref/genome \
--genomeFastaFiles $p_ref/genome/sacCer3.fa

cd $p_trim

# Do alignment
for file in $(ls *_R1*.gz);do
    # Extract the base name without the extension
    base_name="${file%_*_*.*.*}"

    # Construct the input and output file names
    read_1=${base_name}_R1_paired.fastq.gz
    read_2=${base_name}_R2_paired.fastq.gz

    # Align with STAR
    STAR \
    --runThreadN 12 \
    --runMode alignReads \
    --readFilesType Fastx \
    --readFilesCommand zcat \
    --genomeDir $p_ref/genome \
    --sjdbGTFfile $p_ref/annotation/sacCer3.ensGene.gtf \
    --sjdbOverhang 100 \
    --readFilesIn \
    $read_1 \
    $read_2 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMmode Full \
    --outFileNamePrefix $p_align/${base_name}_
done

#move bam file 
cd $p_align
mkdir bam
mv *.bam bam/
