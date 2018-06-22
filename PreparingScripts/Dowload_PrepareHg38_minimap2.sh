
#### HG38 reference ###

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > hg38.fa


### indexing genome for minimap2 (v. 2.10-r768-dirty) ###

minimap2 -d hg38.mmi hg38.fa

