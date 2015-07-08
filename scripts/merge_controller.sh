#!/bin/bash

# initialize
ENV="environment.sh"
if [ -f $ENV ];
then
    echo "Detected environment file!\nRemember that this file defines the reference that should be used in later steps."       
else
    echo "No environment file detected.	The program will now stop!"
    exit 1
fi
source environment.sh

# set variables
pattern=$1
outname=$2
threads=$3
name="$outname"_merged

#########################################

for i in $pattern; do echo $i; done > "$outname"_to_merge

headerfile=`head -n 1 "$outname"_to_merge`

echo "$(date): merging BAM files..."
samtools merge -@ $threads -h $headerfile "$outname"_merged.bam $pattern 
if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during merge"; exit 1; else echo "$(date): merging BAMs done"; fi


