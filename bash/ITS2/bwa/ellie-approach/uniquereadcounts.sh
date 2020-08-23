#!/bin/bash

for file in *trimmed.sam
do
  #echo "$file"
  grep -c XT:A:U $file
done
