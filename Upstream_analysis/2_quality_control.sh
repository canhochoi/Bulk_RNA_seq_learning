#!/usr/bin/bash

p_sra='/home/lqluan/bulk_RNAseq/week7/sacCer2/sra'
p_raw='/home/lqluan/bulk_RNAseq/week7/sacCer2/raw'

for file in $(ls $p_raw/*gz); do
    fastqc $file -o $p_raw/qc_check/
done
