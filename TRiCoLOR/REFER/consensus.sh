cd $2
find . -name "*.unaligned.fa" -printf '%f\n' > files.txt

i=0
while ((i++)); read -r line
do
    $1 consensus -f fasta -t ont $line
    minimap2 -ax map-ont $3 cs.fa.gz > cs.sam
    samtools view -bS cs.sam > cs.bam && rm cs.sam
    samtools sort cs.bam > Fasta$i.consensus.srt.bam && rm cs.bam
    samtools index Fasta$i.consensus.srt.bam

done < files.txt

find . -not -name "*.srt.bam" -not -name "*.srt.bam.bai" -not -name "*.bed" -delete




