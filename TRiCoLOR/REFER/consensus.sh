#!/bin/bash

cd $1

$2 1 $7 $8 $9 $10 $4 > $3".cs.fa" && rm $4
minimap2 -ax $5 $6 $3".cs.fa" > $3".cs.sam" && rm $3".cs.fa"
samtools view -b $3".cs.sam"  > $3".cs.bam" && rm $3".cs.sam"
samtools sort $3".cs.bam" > $3".cs.srt.bam" && rm $3".cs.bam"
samtools index $3".cs.srt.bam"
