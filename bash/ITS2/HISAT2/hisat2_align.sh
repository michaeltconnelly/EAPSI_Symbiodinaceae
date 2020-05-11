#!/bin/bash
#local alignment of ITS2 sequences using HISAT2
EAPSIsamples="Wt2-6b"
for EAPSIsample in $EAPSIsamples
do
/Users/mikeconnelly/computing/programs/hisat2-2.1.0/hisat2 \
-x ./data/indices/HISAT2/symrefs \
-U ./data/reads/${EAPSIsample}_PdamUnmapped.out.mate1.fastq \
-q \
--phred33 \
--score-min L,0.0,-0.4 \
-S ./outputs/symbiont_types/HISAT2/${EAPSIsample}_HISAT2_symrefs_ITS2.sam \
--summary-file ./outputs/symbiont_types/HISAT2/${EAPSIsample}_symrefs_summary.txt
done
