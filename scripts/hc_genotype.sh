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
name="$1"

###########################################

echo "$(date): GATK generating VCF..."
java -Xmx8g -Djava.io.tmpdir="$scratchdir" -jar $GATK_HOME/GenomeAnalysisTK.jar \
-R "$reference" \
-T HaplotypeCaller \
-I "$name" \
-o "$name".gatk.vcf \
--variant_index_type LINEAR \
--variant_index_parameter 128000;
if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during VCF generation GATK"; exit 1; else echo "$(date): VCF generation done"; fi

