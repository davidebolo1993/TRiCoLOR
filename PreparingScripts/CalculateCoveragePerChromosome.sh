#!/bin/bash



File="TestBed.txt" #bed with tandem repeats deletions (HOM/HAP0/HAP1) and fake region in which no alteretions where detected and are randomly distributed

# rescale number of alterations to pass to samtools depth using chromosomes size
# requires bedops

sort-bed $File> elements.bed

for chr in `bedextract --list-chr elements.bed`

do 

bedextract $chr elements.bed > elements.$chr.bed
count=`wc -l elements.$chr.bed | cut -d' ' -f1`
echo -e "$chr\t$count" 

done > counts.txt

sampleSize=300
sum=`cut -f2 counts.txt | perl -nle '$sum += $_ } END { print $sum' -`
awk -v sum=$sum -v sampleSize=$sampleSize '{print $0"\t"($2/sum)"\t"int(sampleSize*$2/sum)}' counts.txt > proportions.txt
awk '{ \
    perChrName=$1; \
    perChrN=$4; \
    cmd="shuf elements."perChrName".bed | head -n "perChrN; \
    system(cmd); \
}' proportions.txt \
| sort-bed - \
> sample.bed

rm elements.chr*
rm counts.txt
rm proportions.txt

# Sample bed file is used to calculate coverage at intervals
# cd /home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled

samtools depth -b /home/davideb/sample.bed ONTReads.bam > Depth.txt 

#coverage per chromosome 

Chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"

OutputDir="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/ChromosomesCoverage"

mkdir $OutputDir

for Index in $Chromosomes

do
	awk -v i="$Index" '$1 == i {print $0}' Depth.txt > $OutputDir"/"$Index.coverage
  
done


# per chromosome coverage results will be analyzed and plotted in python
