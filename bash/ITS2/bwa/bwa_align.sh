#!/bin/bash
bwa aln \
-n 0.005 \
-k 5 \
./data/indices/bwaindex/GeoSymbio \
./data/reads/Wt2-6b_PdamUnmapped.out.mate1.fastq > ./data/reads/Wt2-6b_PdamUnmapped.out.mate1.sai
###
bwa samse ./data/indices/bwaindex/GeoSymbio \
./data/reads/Wt2-6b_PdamUnmapped.out.mate1.sai ./data/reads/Wt2-6b_PdamUnmapped.out.mate1.fastq > ./outputs/ITS2aln/aln-se.sam
