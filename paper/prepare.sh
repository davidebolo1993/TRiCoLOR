#!/bin/bash

##resolving haplotypes of a long-read BAM using publicly available data from Chaisson et al., 2019.
##for in-house long reads sequencing experiment, have a look at longshot (https://github.com/pjedge/longshot)
##if dealing with a haplotype-tagged BAM from whatshap (https://whatshap.readthedocs.io/en/latest/), bamtools should serve the purpose (bamtools split -tag HP)

##requires samtools, bgzip, tabix, bcftools and alfred (https://github.com/tobiasrausch/alfred)
##run with sh prepare.sh in.vcf.gz in.srt.bam in.ref.fa in.alfredpath

##1
##get the phased vcf.gz
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20170323_Strand-seq_phased_FB+GATK_VCFs/HG00733_StrandS_SNVsetFP+GATK_phased.vcf.gz

##2
##get the corresponding BAM (from Pacific Biosciences, for example). They are splitted by chromosome, we need to merge these
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180102_pacbio_blasr_reheader/HG00733*.bam 
##assuming no other BAM in folder
#ls *.bam > bamlist.txt
#samtools merge -@ 20 -m 5G -b bamlist.txt HG00733.srt.bam && samtools index HG00733.srt.bam

##3
##get the reference genome
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
#samtools index GRCh38_full_analysis_set_plus_decoy_hla.fa

##4
##assuming alfred executable is in /opt/alfred/bin/alfred

#sh prepare.sh HG00733_StrandS_SNVsetFP+GATK_phased.vcf.gz HG00733.srt.bam GRCh38_full_analysis_set_plus_decoy_hla.fa /opt/alfred/bin/alfred

vcfgz=${1}
bam=${2}
reference=${3}
alfred=${4}
id=$(echo ${vcfgz} | cut -f1 -d "_")
zcat ${vcfgz} > ${id}.vcf
sed '/##FORMAT=<ID=P2,Number=1,Type=Float,Description="Probability value of allele 2">/a ##SAMPLE=<ID=${id}>' ${id}.vcf > ${id}.fixed.tmp.vcf
sed 's/[[:blank:]]*$//' ${id}.fixed.tmp.vcf > ${id}.fixed.vcf && rm ${id}.vcf && rm ${id}.fixed.tmp.vcf
bgzip ${id}.fixed.vcf && tabix ${id}.fixed.vcf.gz
bcftools view -m2 -M2 -c 1 -C 1 -e 'GT[*] = "mis"' -O b -o ${id}.fixed.bcf ${id}.fixed.vcf.gz && rm ${id}.fixed.vcf.gz*
bcftools index ${id}.fixed.bcf
${alfred} split -r ${reference} -s ${id} -v index ${id}.fixed.bcf -p ${id}.h1.bam -q ${id}.h2.bam ${bam}

#end 
