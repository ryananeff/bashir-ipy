#!/bin/bash

# initialize
ENV="environment.sh"
if [ -f $ENV ];
then
    echo "Detected environment file!"
    echo "Remember that this file defines the reference that should be used in later steps."       
else
    echo "No environment file detected.	The program will now stop!"
    exit 1
fi
source environment.sh

# set variables
input="$1"
threads="$2"
name="$3"

#import variables from the environment file
suffix=$suffix
reference=$reference

###########################################

# initialize
echo "$(date): BLASR aligning file; name: $name..."

blasr $input $reference -nproc $threads -sa $suffix -sam -out $name -clipping subread -bestn 2

#end
if [ $? -ne 0 ]; 
then 
	echo "$(date): exited with non-zero status ($?) during BLASR alignment; name $name"; 
	exit 1; 
else 
	echo "$(date): BLASR alignment done; name $name"; 
fi

# done (yes!)


