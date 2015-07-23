#!/bin/bash

# initialize
ENV="environment.sh"
if [ -f $ENV ];
then
    echo -ne "Detected environment file!\nRemember that this file defines the reference that should be used in later steps.\n"       
else
    echo "No environment file detected.	The program will now stop!"
    exit 1
fi
source environment.sh

# set variables
multivcf=$1 # this should include all samples together in one VCF with individual names in the cols
pedfile=$2  # this needs to have the same names as the sample names in the VCF
outname=$3

###########################################

echo "$(date): GATK generating VCF..."
java -Xmx24g -Djava.io.tmpdir="$scratchdir" -jar $GATK_HOME/GenomeAnalysisTK.jar \
-R "$reference" \
-T PhaseByTransmission \
-V "$multivcf" \
-ped $pedfile \
-o "$outname".phased-trio.vcf ;
if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during VCF phase by trans GATK"; exit 1; else echo "$(date): VCF phase by trans (GATK) done."; fi

