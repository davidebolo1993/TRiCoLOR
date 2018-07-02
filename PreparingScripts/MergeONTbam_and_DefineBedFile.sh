#!/bin/bash

#merge ".srt.bam" files from single .fastq ONT rebasecalled files

#create list of bam file to merge

### ListBam ###

#HG00733_lib01_20171205_FAH36618_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib01_20171205_FAH36664_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib01_20171205_FAH36718_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib02_20171207_FAH36590_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib02_20171207_FAH36591_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib02_20171207_FAH36626_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib02_20171207_FAH36653_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib03_20171207_FAH36732_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib04_20171213_FAH36433_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib04_20171213_FAH37423_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib04_20171213_FAH37574_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib05_20171213_FAH36467_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib05_20171213_FAH36476_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib05_20171213_FAH36782_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib06_20171218_FAH36588_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib06_20171218_FAH36624_DD_guppy_0.5.1.minimap2.srt.bam
#HG00733_lib06_20171218_FAH36752_DD_guppy_0.5.1.minimap2.srt.bam


samtools merge -b ListBam ONTReads.bam
samtools index ONTReads.bam


BedFile="/home/davideb/HG00733.calls_plus_loci.bed"


#Will be modified to extract only deletions
