#!/bin/bash

for i in *.fofn; do 

submitjob 48 -c 32 -m 40 -P acc_HuPac -q premium ~/work/documents/scripts/blasr_align.sh \
$i 32 "$i"_aligned;

done
