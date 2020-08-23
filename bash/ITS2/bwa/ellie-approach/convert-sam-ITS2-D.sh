#!/bin/bash

for file1 in *_trimmed.sai
do
  for file2 in *_trimmed.fasta
  do
    echo "converting $file1"
    bwa samse ITS2-D $file1 $file2 > ${file1:0:6}_trimmed.sam
  done
done
