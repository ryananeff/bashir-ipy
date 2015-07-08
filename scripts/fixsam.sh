#!/bin/bash

infile=$1

cat $infile | sed 's/ RG=/-RG=/g' | samtools view -Shb -o "$infile".fixed.bam -


