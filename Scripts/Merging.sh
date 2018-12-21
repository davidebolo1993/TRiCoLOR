#!/bin/bash

#merging .bam files without exeeding limit of 1024 open files at time

file=$(ls $1'*.srt.bam' | sort -R | tail -1)
samtools view -bH $file > $2.bam
find $1 -name '*srt.bam' -exec samtools view -b {} \; >> $2'.bam'
samtools sort $2'.bam' > $2'.srt.bam'
rm $2'.bam'
samtools index $2'.srt.bam'
