#!/bin/bash

#purpose: quantify non-coral trimmed RNAseq reads against the Cladocopium C1 reference transcriptome (Levin et al. 2016)
#conda activate salmon
salmon index -t ./data/refs/GSE72763_MI_min250_nr.fasta.gz -i ./data/indices/symC1_MI_index

salmon quant -i ./data/indices/symC1_MI_index \
-l SR \
-r ./data/reads/Wt2-6b_PdamUnmapped.out.mate1.fastq \
--validateMappings -o ./outputs/transcripts_quant
