#!/bin/bash

#merging .bam files without exceeding the number of opened files at the same time

cd $2
samtools view -H $1 > $3.sam
find . -type f -name '*srt.bam' -exec samtools view {} \; >> $3'.sam'
htsbox samview -bS $3'.sam' > $3'.bam' && rm $3'.sam'
samtools sort $3'.bam' > $3'.srt.bam'
rm $3'.bam' && samtools index $3'.srt.bam'
find . -type f -not -name '*.merged.srt.bam' -not -name '*merged.srt.bam.bai' -not -name '*.tsv' -delete
