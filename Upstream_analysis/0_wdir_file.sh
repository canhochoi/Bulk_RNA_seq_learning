#!usr/bin/bash

mkdir sacCer2
mkdir -p sacCer2/sra/
mkdir -p sacCer2/ref/
mkdir -p sacCer2/ref/annotation/
mkdir -p sacCer2/ref/genome/
mkdir -p sacCer2/raw/
mkdir -p sacCer2/raw/qc_check
mkdir -p sacCer2/trim/
mkdir -p sacCer2/trim/qc_check
mkdir -p sacCer2/align
mkdir -p sacCer2/align/bam
mkdir -p sacCer2/counts/



cd sacCer2/sra
start=23867632
end=23867637
output_file="SraAccList.csv"

for ((i=start; i<=end; i++))
do
    echo "SRR$i" >> "$output_file"
done
