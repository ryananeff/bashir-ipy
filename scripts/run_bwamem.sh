#!/bin/bash

# initialize
source environment.sh

# set variables
left=$1
right=$2
sample=$3
lane=$4
study=$5
name="$sample"_"$lane"

###########################################

# run new bwamem alignment, use two threads and 8g of RAM (we're being hopeful here...)
#echo "$(date): starting bwamem alignment...";
bwa mem -t 12 -w 100 -k 20 -M -R '@RG\tID:1\tSM:$3\tLB:$3\tPL:ILLUMINA\tPU:$3' "$reference" "$left" "$right" > "$name".sam;
#if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during bwa-mem"; exit 1; else echo "$(date): bwamem alignment done."; fi

# create bam file and sort it before writing
#echo "$(date): starting conversion to BAM...";
#samtools view -bS "$name".sam -o "$name".bam;
#if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during samtools view"; exit 1; else echo "$(date): BAM conversion done."; fi

#echo "$(date): sorting BAM...";
#samtools sort "$name".bam "$name".sort;
#if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during samtools sort"; exit 1; else echo "$(date): sorting bam done."; fi

# if all went well, remove the old files (we should only get here if there was zero exit statuses
#if [ -s "$name".sort.bam ]
#then echo "$(date): sorted BAM file found with non-zero size; continuing..."
#else echo "$(date): unexpected error: sorted BAM file is missing; exiting"; exit 2; fi
#echo "$(date): removing temporary files..."
#rm "$name".sam "$name".bam
#if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during temporary file removal"; exit 1; else echo "$(date): removed temporary SAM and BAM."; fi

# index bam again
#echo "$(date): indexing sorted BAM..."
#samtools index "$name".sort.bam;
#if [ $? -ne 0 ]; then echo "$(date): exited with non-zero status ($?) during samtools index sorted"; exit 1; else echo "$(date): indexing sorted bam done."; fi

# now we need to merge by lane (can't be handled inside this program)

