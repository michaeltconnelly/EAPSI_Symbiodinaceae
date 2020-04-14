#!/bin/bash
#local alignment of ITS2 sequences using HISAT2
EAPSIsamples="Wt1-6a Wt1-6b Wt1-6c Wt2-6a Wt2-6b Wt2-6c Hw1-6a Hw1-6b Hw1-6c Hw2-6b Hw2-6c"
for EAPSIsample in $EAPSIsamples
do
/Users/mikeconnelly/computing/programs/hisat2-2.1.0/hisat2 \
-x ./data/indices/HISAT2index/SymITS2 \
-U ./data/reads/${EAPSIsample}_PdamUnmapped.out.mate1.fastq \
-q \
--phred33 \
--score-min L,0.0,-0.4 \
-S ./outputs/HISAT2align_ITS2/${EAPSIsample}_HISAT2_ITS2.sam \
--summary-file ./outputs/HISAT2align_ITS2/${EAPSIsample}_summary.txt
done
