#!/bin/bash

echo "Checking TRiCoLOR installation"
echo ""
echo ""
TRiCoLOR -h
echo ""
echo ""
echo "Testing TRiCoLOR SENSoR"
echo ""
echo ""
TRiCoLOR SENSoR -h
TRiCoLOR SENSoR -bam son/sim.srt.bam -o sensor_son
echo ""
echo ""
echo "Testing TRiCoLOR REFER"
echo ""
echo ""
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
TRiCoLOR REFER -h
TRiCoLOR REFER -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bam son/sim.srt.bam -bed sensor_son/TRiCoLOR.srt.bed.gz -o refer_son -s 20 --samplename SON
echo ""
echo ""
echo "Testing TRiCoLOR ApP"
echo ""
echo ""
TRiCoLOR ApP -h
TRiCoLOR ApP -g GRCh38_full_analysis_set_plus_decoy_hla.fa -bam refer_son/haplotype1/TRiCoLOR.srt.bam refer_son/haplotype2/TRiCoLOR.srt.bam -o app_son -gb refer_son/reference/TRiCoLOR.srt.bed.gz -h1b refer_son/haplotype1/TRiCoLOR.srt.bed.gz -h2b refer_son/haplotype2/TRiCoLOR.srt.bed.gz chr20:17553795-17553825
echo ""
echo ""
echo "Testing TRiCoLOR SAGE"
echo ""
echo ""
TRiCoLOR SAGE -h
TRiCoLOR SAGE -vcf refer_son/TRiCoLOR.srt.vcf.gz -bam parent1/sim.srt.bam parent2/sim.srt.bam -o sage_trio --samplename PARENT1 PARENT2 --mendel