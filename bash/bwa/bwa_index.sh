#!/bin/bash
#purpose: create BWA index from ITS2 reference

#BSUB -J bwa_index
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o bwa_idx%J.out
#BSUB -e bwa_idx%J.err
#BSUB -n 8
#BSUB -W 6:00
#BSUB -u mconnelly@rsmas.miami.edu
#BSUB -N

mkdir ./data/indices/bwa_index/

module load bwa/0.7.4
bwa index \
-p ./data/bwa_index/ \
-a is \
./data/GeoSymbio/Symbiodatabaceae.fasta
