#!/bin/bash
EAPSIsamples="Wt1-6a Wt1-6b Wt1-6c Wt2-6a Wt2-6b Wt2-6c Hw1-6a Hw1-6b Hw1-6c Hw2-6b Hw2-6c"
for EAPSIsample in $EAPSIsamples
do
samtools view -bS ./outputs/HISAT2align_ITS2/${EAPSIsample}_HISAT2_ITS2.sam > ${EAPSIsample}_HISAT2_ITS2.bam
samtools sort ${EAPSIsample}_HISAT2_ITS2.bam -o ./outputs/HISAT2align_ITS2/${EAPSIsample}_HISAT2_ITS2.sorted.bam
samtools index -b ./outputs/HISAT2align_ITS2/${EAPSIsample}_HISAT2_ITS2.sorted.bam ./outputs/HISAT2align_ITS2/${EAPSIsample}_HISAT2_ITS2.sorted.bai
done


samtools view -F 4 in.sam > mapped.sam
