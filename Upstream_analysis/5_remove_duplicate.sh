#!usr/bin/bash
# Create path to your picard tool
picard_jar='/home/lqluan/bulk_RNAseq/week7/tools/picard/build/libs/picard.jar'
p_align='/home/lqluan/bulk_RNAseq/week7/sacCer2/align'

cd $p_align/bam/

# Using for loop to remove all duplicate in bam file
for file in $(ls *.bam);do
    base_name="${file%Aligned*}"
    java -jar $picard_jar MarkDuplicates \
    --INPUT ${base_name}Aligned.sortedByCoord.out.bam \
    --OUTPUT ${base_name}_rmdup.bam \
    --METRICS_FILE ${base_name}_rmdup.metrics2 \
    --REMOVE_DUPLICATES true
done
