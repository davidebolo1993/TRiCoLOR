#!/bin/bash

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > hg38.fa

### indexing for Minimap2

minimap2 -d hg38.mmi hg38.fa
