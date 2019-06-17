#!/bin/bash

cd $1

if [ ! -z `find . -type f -name "*.srt.bam" -not -name "chr*" -print -quit` ]; then

	samtools view -H $2 > $3".merged.sam"
	find . -type f -name "*.srt.bam" -not -name "chr*" -exec samtools view {} \; >> $3".merged.sam"
	samtools view -b $3".merged.sam" > $3".merged.bam" && rm $3".merged.sam"
	samtools sort -@ $4 $3".merged.bam" > $3".merged.srt.bam"
	rm $3".merged.bam" && samtools index $3".merged.srt.bam"
	find . -type f -not -name "*.merged.srt.bam" -not -name "*.merged.srt.bam.bai" -not -name "*.bed" -delete

fi
