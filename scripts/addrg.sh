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
name=$1
lane=$2
sample=$3
study=$4

###########################################

# add back read groups
echo "$(date): adding back read groups..."
java -Xmx8g -Djava.io.tmpdir="$scratchdir" -jar $picardhome/AddOrReplaceReadGroups.jar \
I="$name" \
O="$name".rg.bam \
RGID="$lane" \
RGLB="$sample" \
RGPL=Illumina \
RGSM="$sample" \
RGPU="$study" \
VALIDATION_STRINGENCY=LENIENT;
if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during read group adding"; exit 1; else echo "$(date): adding readgroups done"; fi

# index bam again
echo "$(date): indexing rg BAM..."
samtools index "$name".rg.bam
if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during samtools index rg"; exit 1; else echo "$(date): indexing rg bams done."; fi

# now we need to merge by lane (can't be handled inside this program)

