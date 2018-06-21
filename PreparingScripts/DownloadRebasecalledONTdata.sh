#!/bin/bash

while true;do
wget -T 15 -r -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20180502_ONT_rebasecalled/ && break
done

