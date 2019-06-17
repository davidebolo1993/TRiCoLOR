cd $1

$2 consensus -f fasta -t $3 -a $4".al.fa.gz" -c $4".cs.fa.gz" $5 && rm $5 && rm $5".fai"
minimap2 -ax $6 $7 $4".cs.fa.gz" > $4".cs.sam" && rm $4".al.fa.gz" && rm $4".cs.fa.gz"
samtools view -b $4".cs.sam" > $4".cs.bam" && rm $4".cs.sam"
samtools sort $4".cs.bam" > $4".cs.srt.bam" && rm $4".cs.bam"
samtools index $4".cs.srt.bam"
