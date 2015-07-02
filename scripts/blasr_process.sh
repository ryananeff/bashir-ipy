#!/bin/bash

for i in *.fofn; do 

submitjob 12 -c 12 -m 31 -P acc_HuPac -q premium ~/work/documents/scripts/blasr_align.sh \
$i 12 "$i"_aligned;

done
