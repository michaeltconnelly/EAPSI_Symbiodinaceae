#!/bin/bash

for file in *_trimmed.fasta
do
  echo "converting $file"
  bwa aln -n 0.005 -k 5 ITS2-D $file > ${file:0:6}_trimmed.sai
done
