#!/bin/bash

# initialize
ENV="environment.sh"
if [ -f $ENV ];
then
    echo -ne "Detected environment file!\nRemember that this file defines the reference that should be used in later steps."       
else
    echo "No environment file detected.	The program will now stop!"
    exit 1
fi
source environment.sh

# set variables
name="$1"
infile="$2"

###########################################

echo "$(date): freebayes generating VCF..."
freebayes -f $reference \
--genotype-qualities \
--report-genotype-likelihood-max \
--standard-filters \
--report-all-haplotype-alleles \
-b $infile \
-v "$name".freebayes.vcf;
if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during VCF generation freebayes"; exit 1; else echo "$(date): Freebayes VCF generation done"; fi

