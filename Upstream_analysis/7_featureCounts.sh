#!/usr/bin/bash
p_ref='/home/lqluan/bulk_RNAseq/week7/sacCer2/ref'
p_align='/home/lqluan/bulk_RNAseq/week7/sacCer2/align'
p_trim='/home/lqluan/bulk_RNAseq/week7/sacCer2/trim'
p_bam='/home/lqluan/bulk_RNAseq/week7/sacCer2/align/bam'
p_annotation='/home/lqluan/bulk_RNAseq/week7/sacCer2/ref/annotation'
p_counts='/home/lqluan/bulk_RNAseq/week7/sacCer2/counts'

featureCounts \
-p \
-a $p_annotation/sacCer3.ensGene.gtf \
-o $p_counts/all.featureCounts.txt \
$p_bam/WT_C_1_rmdup.bam \
$p_bam/WT_C_2_rmdup.bam \
$p_bam/WT_C_3_rmdup.bam \
$p_bam/WT_E_1_rmdup.bam \
$p_bam/WT_E_2_rmdup.bam \
$p_bam/WT_E_3_rmdup.bam


