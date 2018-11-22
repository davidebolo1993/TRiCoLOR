#!/bin/bash

############################################ Download ONT Data (HG00733) ############################################

while true;do
wget -T 15 -r -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/ && break
done


############################################ Download Hg38 Reference ############################################

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > hg38.fa


############################################ Minimap2 Index for Hg38 Reference ############################################

#minimap2 2.11


Minimap2="/home/davideb/anaconda2/pkgs/minimap2-2.11-ha92aebf_0/bin/minimap2" 



$Minimap2 -d hg38.mmi hg38.fa



############################################ Calculate Minimap2 Performances ############################################


### Time Performances ###

find . -maxdepth 1 -name "*gz"  | awk '{
	res=$1
	split(res,res2,"./")
	print res2[2]
}'> FileNames.txt


find . -maxdepth 1 -name "*.gz" | while read -r file; do zcat -f "$file" | wc -l ; done > FileLengths.txt #has to be divided by 4 as each sequence is described by 4 lines in .fastq file

# Without trimming on length/quality

# samtools 1.7
# htsbox r340


MyFiles="HG00733_lib01_20171205_FAH36618_DD_guppy_0.5.1.fq.gz HG00733_lib01_20171205_FAH36664_DD_guppy_0.5.1.fq.gz HG00733_lib01_20171205_FAH36718_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36590_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36591_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36626_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36653_DD_guppy_0.5.1.fq.gz HG00733_lib03_20171207_FAH36732_DD_guppy_0.5.1.fq.gz HG00733_lib04_20171213_FAH36433_DD_guppy_0.5.1.fq.gz HG00733_lib04_20171213_FAH37423_DD_guppy_0.5.1.fq.gz HG00733_lib04_20171213_FAH37574_DD_guppy_0.5.1.fq.gz HG00733_lib05_20171213_FAH36467_DD_guppy_0.5.1.fq.gz HG00733_lib05_20171213_FAH36476_DD_guppy_0.5.1.fq.gz HG00733_lib05_20171213_FAH36782_DD_guppy_0.5.1.fq.gz HG00733_lib06_20171218_FAH36588_DD_guppy_0.5.1.fq.gz HG00733_lib06_20171218_FAH36624_DD_guppy_0.5.1.fq.gz HG00733_lib06_20171218_FAH36752_DD_guppy_0.5.1.fq.gz"

Hg38mmRef="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/hg38Reference/hg38.mmi"

FileTime="time.txt"

Time="/usr/bin/time"

Minimap2="/home/davideb/anaconda2/pkgs/minimap2-2.11-ha92aebf_0/bin/minimap2" 


for Index in $MyFiles

do

zcat $Index | awk '{if(NR%4==2) print length($1)}' > $Index"_Length.txt" && #length of each sequence
cat $Index"_Length.txt" | awk 'BEGIN{FS=OFS="\t"; }{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1< min) {min=$1}; total+=$1; count+=1} END {print "mean =", total/count, "\nminimum =", min, "\nmaximum =", max}'> $Index"_MeanMinMaxLength.txt" &&
cat $Index"_MeanMinMaxLength.txt" | grep "mean" | awk '{print $3}' > $Index"_MeanLength.txt" &&
cat $Index"_MeanMinMaxLength.txt" | grep "minimum" | awk '{print $3}' > $Index"_MinLength.txt" &&
cat $Index"_MeanMinMaxLength.txt" | grep "maximum" | awk '{print $3}' > $Index"_MaxLength.txt" &&
rm $Index"_Length.txt" $Index"_MeanMinMaxLength.txt" &&
$Time -p -ao $FileTime $Minimap2 -ax map-ont -t 10 $Hg38mmRef $Index > $Index".minimap2.sam" &&
htsbox samview -bS $Index.minimap2.sam > $Index".minimap2.bam" &&
rm $Index.minimap2.sam &&
samtools sort -m 5G -@ 10 $Index.minimap2.bam > $Index".minimap2.srt.bam" &&
rm $Index.minimap2.bam &&
samtools index $Index.minimap2.srt.bam &&
mv $Index.minimap2.srt.bam "$(echo "$Index.minimap2.srt.bam" | sed -e 's/.fq.gz//')" &&
mv $Index.minimap2.srt.bam.bai "$(echo "$Index.minimap2.srt.bam.bai" | sed -e 's/.fq.gz//')"

done


cat time.txt | grep "real" | awk '{print $2}' > RealTime.txt
cat *_MeanLength.txt > MeanLength.txt
rm *_MeanLength.txt
cat *_MinLength.txt > MinLength.txt
rm *_MinLength.txt
cat *_MaxLength.txt > MaxLength.txt
rm *_MaxLength.txt


paste -d "\t" FileNames.txt FileLengths.txt RealTime.txt MeanLength.txt MinLength.txt MaxLength.txt > Informations.txt 

awk 'BEGIN{FS=OFS="\t" }{$2 = $2 /4}1{$3 = $3 /60}1' Informations.txt > Minimap2_PyInfo.txt #true length of fastq and convert seconds to minutes

rm Informations.txt FileNames.txt FileLengths.txt RealTime.txt MeanLength.txt MinLength.txt MaxLength.txt time.txt

mkdir Minimap2_Time_Performances 

cp Minimap2_PyInfo.txt *.srt.bam* Minimap2_Time_Performances/ &&

rm Minimap2_PyInfo.txt *.srt.bam*



### Accuracy in Alignment ###

# pbsim v1.0.3 


mkdir Minimap2_Accuracy_Performances

cd Minimap2_Accuracy_Performances

#sample from hg38 

Hg38Dir="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/hg38Reference"


samtools faidx $Hg38Dir"/hg38.fa" chr1:1000001-40000000 > $Hg38Dir"/Hg38selected.fa" ### 40Mb ##

ReferenceTosample="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/hg38Reference/Hg38selected.fa"

#model qc taken from github page of Pbsim

Modelqc="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2_Accuracy_Performances/model_qc_clr"

#Default values
#Model-based simulation
# CLR == Continuous Long Reads
#--depth 20.0
#--length-max 25000
#--length-min 100
#--length-mean 3000.0
#--accuracy-min 0.75
#--accuracy-max 1.00
#--accuracy-mean 0.78

#List Of abbreviations

#SR = ShorterReads
#MR = MediumReads
#LR = LongerReads
#LA = LowAccuracy
#MA = MediumAccuracy
#HA = HighAccuracy

# ratio 10:60:30 for substitution:insertion:deletion

#random dataset

pbsim --prefix SR_LA --data-type CLR --length-min 100 --length-mean 3000 --length-max 25000 --accuracy-min 0.80 --accuracy-mean 0.85 --accuracy-max 0.90 --model_qc $Modelqc $ReferenceTosample &&
pbsim --prefix MR_LA --data-type CLR --length-min 500 --length-mean 8000 --length-max 35000 --accuracy-min 0.80 --accuracy-mean 0.85 --accuracy-max 0.90 --model_qc $Modelqc $ReferenceTosample &&
pbsim --prefix LR_LA --data-type CLR --length-min 1000 --length-mean 13000 --length-max 45000 --accuracy-min 0.80 --accuracy-mean 0.85 --accuracy-max 0.90 --model_qc $Modelqc $ReferenceTosample &&


pbsim --prefix SR_MA --data-type CLR --length-min 100 --length-mean 3000 --length-max 25000 --accuracy-min 0.85 --accuracy-mean 0.90 --accuracy-max 0.95 --model_qc $Modelqc $ReferenceTosample &&
pbsim --prefix MR_MA --data-type CLR --length-min 500 --length-mean 8000 --length-max 35000 --accuracy-min 0.85 --accuracy-mean 0.90 --accuracy-max 0.95 --model_qc $Modelqc $ReferenceTosample &&
pbsim --prefix LR_MA --data-type CLR --length-min 1000 --length-mean 13000 --length-max 45000 --accuracy-min 0.85 --accuracy-mean 0.90 --accuracy-max 0.95 --model_qc $Modelqc $ReferenceTosample &&


pbsim --prefix SR_HA --data-type CLR --length-min 100 --length-mean 3000 --length-max 25000 --accuracy-min 0.90 --accuracy-mean 0.95 --accuracy-max 0.99 --model_qc $Modelqc $ReferenceTosample &&
pbsim --prefix MR_HA --data-type CLR --length-min 500 --length-mean 8000 --length-max 35000 --accuracy-min 0.90 --accuracy-mean 0.95 --accuracy-max 0.99 --model_qc $Modelqc $ReferenceTosample &&
pbsim --prefix LR_HA --data-type CLR --length-min 1000 --length-mean 13000 --length-max 45000 --accuracy-min 0.90 --accuracy-mean 0.95 --accuracy-max 0.99 --model_qc $Modelqc $ReferenceTosample &&


# dataset with min/mean/max similar to ONT rebasecalled data

#accuracy 0.8-0.85-0.9

pbsim --prefix ONT --data-type CLR --length-min 100 --length-mean 10000 --length-max 150000 --accuracy-min 0.80 --accuracy-mean 0.90 --accuracy-max 0.99 --model_qc $Modelqc $ReferenceTosample &&

# dataset with high high accuracy, for ultra-long reads

#LLR = LongLongReads
#HHA = HighHighAccurcy 

pbsim --prefix LLR_HHA --data-type CLR --length-min 50000 --length-mean 100000 --length-max 200000 --accuracy-min 0.97 --accuracy-mean 0.98 --accuracy-max 0.99 --model_qc $Modelqc $ReferenceTosample



MyFiles="SR_LA MR_LA LR_LA SR_MA MR_MA LR_MA SR_HA MR_HA LR_HA ONT LLR_HHA"

Hg38mmRef="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/hg38Reference/hg38.mmi"

Minimap2="/home/davideb/anaconda2/pkgs/minimap2-2.11-ha92aebf_0/bin/minimap2" 



for Index in $MyFiles

do 

$Minimap2 -ax map-ont -t 10 $Hg38mmRef $Index"_0001.fastq" > $Index".minimap2.sam" &&
htsbox samview -bS $Index.minimap2.sam > $Index".minimap2.bam" &&
rm $Index.minimap2.sam &&
samtools sort -m 5G -@ 10 $Index.minimap2.bam > $Index".minimap2.srt.bam" &&
rm $Index.minimap2.bam &&
samtools index $Index.minimap2.srt.bam

done


# See .py codes for plotting in Python



############################################ Calculate NGMLR Performances ############################################

#ngmlr 0.2.7

#reference index already created


MyFiles="HG00733_lib01_20171205_FAH36618_DD_guppy_0.5.1.fq.gz HG00733_lib01_20171205_FAH36664_DD_guppy_0.5.1.fq.gz HG00733_lib01_20171205_FAH36718_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36590_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36591_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36626_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36653_DD_guppy_0.5.1.fq.gz HG00733_lib03_20171207_FAH36732_DD_guppy_0.5.1.fq.gz HG00733_lib04_20171213_FAH36433_DD_guppy_0.5.1.fq.gz HG00733_lib04_20171213_FAH37423_DD_guppy_0.5.1.fq.gz HG00733_lib04_20171213_FAH37574_DD_guppy_0.5.1.fq.gz HG00733_lib05_20171213_FAH36467_DD_guppy_0.5.1.fq.gz HG00733_lib05_20171213_FAH36476_DD_guppy_0.5.1.fq.gz HG00733_lib05_20171213_FAH36782_DD_guppy_0.5.1.fq.gz HG00733_lib06_20171218_FAH36588_DD_guppy_0.5.1.fq.gz HG00733_lib06_20171218_FAH36624_DD_guppy_0.5.1.fq.gz HG00733_lib06_20171218_FAH36752_DD_guppy_0.5.1.fq.gz"

Hg38Ref="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/hg38Reference/hg38.fa"

FileTime="time.txt"

Time="/usr/bin/time"

NGMLR="/home/davideb/anaconda2/pkgs/ngmlr-0.2.7-he941832_0/bin/ngmlr"


for Index in $MyFiles

do


$Time -p -ao $FileTime $NGMLR -t 10 -r $Hg38Ref -q $Index -o $Index".ngmlr.sam" -x ont &&
htsbox samview -bS $Index.ngmlr.sam > $Index".ngmlr.bam" &&
rm $Index.ngmlr.sam &&
samtools sort -m 5G -@ 10 $Index.ngmlr.bam > $Index".ngmlr.srt.bam" &&
rm $Index.ngmlr.bam &&
samtools index $Index.ngmlr.srt.bam &&
mv $Index.ngmlr.srt.bam "$(echo "$Index.ngmlr.srt.bam" | sed -e 's/.fq.gz//')" &&
mv $Index.ngmlr.srt.bam.bai "$(echo "$Index.ngmlr.srt.bam.bai" | sed -e 's/.fq.gz//')"

done

cat time.txt | grep "real" | awk '{print $2}' > RealTime.txt

#Retrive info using the table generated before (see Minimap2 performances)

cp "/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2_Time_Performances/Minimap2_PyInfo.txt" .

awk 'BEGIN{FS=OFS="\t" }FNR==NR{a[NR]=$1;next}{$3=a[FNR]}1' RealTime.txt Minimap2_PyInfo.txt > Informations.txt

awk 'BEGIN{FS=OFS="\t" }{$3 = $3 /60}1' Informations.txt > Ngmlr_PyInfo.txt #Convert seconds to minutes

rm Informations.txt Minimap2_PyInfo.txt time.txt RealTime.txt 

mkdir Ngmlr_Time_Performances 

cp Ngmlr_PyInfo.txt *.srt.bam* Ngmlr_Time_Performances/ &&

rm Ngmlr_PyInfo.txt *.srt.bam*



### Accuracy in Alignment ###

# pbsim v1.0.3 


mkdir Ngmlr_Accuracy_Performances

cd Ngmlr_Accuracy_Performances

#retrive Pbsim simulations from previous generated samples (see Minimap2 performances)

DirFastq="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2_Accuracy_Performances/"

MyFiles="SR_LA MR_LA LR_LA SR_MA MR_MA LR_MA SR_HA MR_HA LR_HA ONT LLR_HHA"

Hg38Ref="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/hg38Reference/hg38.fa"

NGMLR="/home/davideb/anaconda2/pkgs/ngmlr-0.2.7-he941832_0/bin/ngmlr"



for Index in $MyFiles

do 

$NGMLR -t 10 -r $Hg38Ref -q $DirFastq$Index"_0001.fastq" -o $Index".ngmlr.sam" -x ont &&
htsbox samview -bS $Index.ngmlr.sam > $Index".ngmlr.bam" &&
rm $Index.ngmlr.sam &&
samtools sort -m 5G -@ 10 $Index.ngmlr.bam > $Index".ngmlr.srt.bam" &&
rm $Index.ngmlr.bam &&
samtools index $Index.ngmlr.srt.bam

done


# See .py codes for plotting in Python


############################################ Calculate BWA performances ############################################



#bwa 0.7.17-r1188

#bwa genome indexing


Hg38Dir="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/hg38Reference/"

cd $Hg38Dir

bwa index -a bwtsw hg38.fa

#back to working directory


Hg38Ref="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/hg38Reference/hg38.fa"

MyFiles="HG00733_lib01_20171205_FAH36618_DD_guppy_0.5.1.fq.gz HG00733_lib01_20171205_FAH36664_DD_guppy_0.5.1.fq.gz HG00733_lib01_20171205_FAH36718_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36590_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36591_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36626_DD_guppy_0.5.1.fq.gz HG00733_lib02_20171207_FAH36653_DD_guppy_0.5.1.fq.gz HG00733_lib03_20171207_FAH36732_DD_guppy_0.5.1.fq.gz HG00733_lib04_20171213_FAH36433_DD_guppy_0.5.1.fq.gz HG00733_lib04_20171213_FAH37423_DD_guppy_0.5.1.fq.gz HG00733_lib04_20171213_FAH37574_DD_guppy_0.5.1.fq.gz HG00733_lib05_20171213_FAH36467_DD_guppy_0.5.1.fq.gz HG00733_lib05_20171213_FAH36476_DD_guppy_0.5.1.fq.gz HG00733_lib05_20171213_FAH36782_DD_guppy_0.5.1.fq.gz HG00733_lib06_20171218_FAH36588_DD_guppy_0.5.1.fq.gz HG00733_lib06_20171218_FAH36624_DD_guppy_0.5.1.fq.gz HG00733_lib06_20171218_FAH36752_DD_guppy_0.5.1.fq.gz"

FileTime="time.txt"

Time="/usr/bin/time"

for Index in $MyFiles

do

$Time -p -ao $FileTime bwa mem -t 10 -x ont2d $Hg38Ref $Index > $Index".bwa.sam" &&
htsbox samview -bS $Index.bwa.sam > $Index".bwa.bam" &&
rm $Index.bwa.sam &&
samtools sort -m 5G -@ 10 $Index.bwa.bam > $Index".bwa.srt.bam" &&
rm $Index.bwa.bam &&
samtools index $Index.bwa.srt.bam &&
mv $Index.bwa.srt.bam "$(echo "$Index.bwa.srt.bam" | sed -e 's/.fq.gz//')" &&
mv $Index.bwa.srt.bam.bai "$(echo "$Index.bwa.srt.bam.bai" | sed -e 's/.fq.gz//')"

done

cat time.txt | grep "real" | awk '{print $2}' > RealTime.txt

#Retrive info using the table generated before (see Minimap2 performances)

cp "/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2_Time_Performances/Minimap2_PyInfo.txt" .

awk 'BEGIN{FS=OFS="\t" }FNR==NR{a[NR]=$1;next}{$3=a[FNR]}1' RealTime.txt Minimap2_PyInfo.txt > Informations.txt

awk 'BEGIN{FS=OFS="\t" }{$3 = $3 /60}1' Informations.txt > BWA_PyInfo.txt #Convert seconds to minutes

rm Informations.txt Minimap2_PyInfo.txt time.txt RealTime.txt 

mkdir BWA_Time_Performances 

cp BWA_PyInfo.txt *.srt.bam* BWA_Time_Performances/ &&

rm BWA_PyInfo.txt *.srt.bam*


### Accuracy in Alignment ###

# pbsim v1.0.3 


mkdir BWA_Accuracy_Performances

cd BWA_Accuracy_Performances


#retrive Pbsim simulations from previous generated samples (see Minimap2 performances)

DirFastq="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2_Accuracy_Performances/"

MyFiles="SR_LA MR_LA LR_LA SR_MA MR_MA LR_MA SR_HA MR_HA LR_HA ONT LLR_HHA"

Hg38Ref="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/hg38Reference/hg38.fa"


for Index in $MyFiles

do 

bwa mem -t 10 -x ont2d $Hg38Ref $DirFastq$Index"_0001.fastq" > $Index".bwa.sam" &&
htsbox samview -bS $Index.bwa.sam > $Index".bwa.bam" &&
rm $Index.bwa.sam &&
samtools sort -m 5G -@ 10 $Index.bwa.bam > $Index".bwa.srt.bam" &&
rm $Index.bwa.bam &&
samtools index $Index.bwa.srt.bam

done


