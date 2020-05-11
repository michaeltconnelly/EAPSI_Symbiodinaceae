#!/bin/bash
#local alignment using HISAT2
/Users/mikeconnelly/computing/programs/hisat2-2.1.0/hisat2-build  \
 -f \
./data/symtypes/Symbiodatabaceae.fasta \
./data/indices/HISAT2/Symbiodatabaceae

/Users/mikeconnelly/computing/programs/hisat2-2.1.0/hisat2-build  \
 -f \
./data/symtypes/GeoSymbio_ITS2_LocalDatabase.fasta \
./data/indices/HISAT2/GeoSymbio

/Users/mikeconnelly/computing/programs/hisat2-2.1.0/hisat2-build  \
 -f \
./data/symtypes/symrefs.fasta \
./data/indices/HISAT2/symrefs
