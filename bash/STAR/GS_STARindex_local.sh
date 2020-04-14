#!/bin/bash
#purpose: create STAR index from ITS2 reference

#BSUB -J star_index
#BSUB -q general
#BSUB -P transcriptomics
#BSUB -o star_idx%J.out
#BSUB -e star_idx%J.err
#BSUB -n 8
#BSUB -W 6:00
#BSUB -u mconnelly@rsmas.miami.edu

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
prodir="/scratch/projects/transcriptomics/mikeconnelly/projects/EAPSI_symITS2"

/Users/mikeconnelly/computing/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--runMode genomeGenerate \
--runThreadN 8 \
--genomeDir /Users/mikeconnelly/computing/sequences/genomes/symbiodinium/GeoSymbio/STARindex \
--genomeFastaFiles /Users/mikeconnelly/computing/sequences/genomes/symbiodinium/GeoSymbio/Symbiodatabaceae.fasta \
--genomeSAindexNbases 6
