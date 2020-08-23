#!/bin/bash
#purpose: create BWA index from ITS2 reference
bwa index \
-p ./data/indices/bwaindex/Symbiodatabaceae \
-a is \
./data/GeoSymbio/Symbiodatabaceae.fasta

bwa index \
-p ./data/indices/bwaindex/GeoSymbio \
-a is \
./data/GeoSymbio/GeoSymbio_ITS2_LocalDatabase.fasta
