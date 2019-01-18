#!/bin/bash

#merging .bam files without exeeding limit of 1024 open files at time

cd $2

samtools view -H $1 > $3'.merged.sam'
find . -type f -name '*.srt.bam' -not -name 'chr*' -exec samtools view -b {} \; >> $3'.merged.sam'
htsbox samview -bS $3'.merged.sam' > $3'.merged.bam' && rm $3'.merged.sam'
samtools sort $3'.merged.bam' > $3'.merged.srt.bam'
rm $3'.merged.bam' && samtools index $3'.merged.srt.bam'
find . -type f -not -name '*.merged.srt.bam' -not -name '*.merged.srt.bam.bai' -not -name '*.tsv' -delete
