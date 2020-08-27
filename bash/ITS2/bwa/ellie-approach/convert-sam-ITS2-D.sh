#!/bin/bash

for file in *_trimmed.sai
do
   echo "converting $file"
   bwa samse ITS2-D $file ${file:0:6}_trimmed.fasta > ${file:0:6}_trimmed.sam
done
