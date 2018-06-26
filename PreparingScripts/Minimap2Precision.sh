#!/usr/bin/bash


###### Minimap2 precision #####

### using  pbsim v1.0.3 ###


mkdir Minimap2Mapping
cd Minimap2Mapping

#sample from hg38 

samtools faidx /home/davideb/nanopore/hg38Reference/hg38.fa chr1:1000001-40000000 > /home/davideb/nanopore/hg38Reference/Hg38selected.fa ### 40Mb ##

ReferenceTosample="/home/davideb/nanopore/hg38Reference/Hg38selected.fa"
Modelqc="/home/davideb/STRCallerProject/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/Minimap2Mapping/model_qc_clr"


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

## random dataset

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
#

pbsim --prefix ONT --data-type CLR --length-min 100 --length-mean 10000 --length-max 150000 --accuracy-min 0.80 --accuracy-mean 0.90 --accuracy-max 0.99 --model_qc $Modelqc $ReferenceTosample &&


# dataset with high high accuracy, for ultra-long reads

#LLR = LongLongReads
#HHA = HighHighAccurcy 

pbsim --prefix LLR_HHA --data-type CLR --length-min 50000 --length-mean 100000 --length-max 200000 --accuracy-min 0.97 --accuracy-mean 0.98 --accuracy-max 0.99 --model_qc $Modelqc $ReferenceTosample



MyFiles="SR_LA MR_LA LR_LA SR_MA MR_MA LR_MA SR_HA MR_HA LR_HA ONT LLR_HHA"

Hg38mmRef="/home/davideb/nanopore/hg38Reference/hg38.mmi"

for Index in $MyFiles

do 

minimap2 -ax map-ont -t 10 $Hg38mmRef $Index"_0001.fastq" > $Index".minimap2.sam" &&
htsbox samview -bS $Index.minimap2.sam > $Index".minimap2.bam" &&
rm $Index.minimap2.sam &&
samtools sort -m 5G -@ 10 $Index.minimap2.bam > $Index".minimap2.srt.bam" &&
rm $Index.minimap2.bam &&
samtools index $Index.minimap2.srt.bam

done

## will be plotted in python

