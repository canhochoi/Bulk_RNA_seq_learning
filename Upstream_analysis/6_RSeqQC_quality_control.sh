#!/usr/bin/bash
p_ref='/home/lqluan/bulk_RNAseq/week7/sacCer2/ref'
p_align='/home/lqluan/bulk_RNAseq/week7/sacCer2/align'
p_trim='/home/lqluan/bulk_RNAseq/week7/sacCer2/trim'
p_bam='/home/lqluan/bulk_RNAseq/week7/sacCer2/align/bam'

#index the bam
samtools index WT_C_1_rmdup.bam

#gene body coverage
geneBody_coverage.py \
    -r $p_ref/annotation/sacCer3.ens.bed \
    -i $p_bam/WT_C_1_rmdup.bam  \
    -o gene_coverage

#junction annotation
junction_annotation.py \
    -r $p_ref/annotation/sacCer3.ens.bed \
    -i $p_bam/WT_C_1_rmdup.bam \
    -o WT_C_1

#junction saturation
junction_saturation.py \
    -r $p_ref/annotation/sacCer3.ens.bed \
    -i $p_bam/WT_C_1_rmdup.bam \
    -o WT_C_1

#read distribution
read_distribution.py \
    -i $p_bam/WT_C_1_rmdup.bam \
    -r $p_ref/annotation/sacCer3.ens.bed	
