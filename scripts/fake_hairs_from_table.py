#!/usr/bin/env python

import pysam
import HTSeq as ht
import math
import os, sys, getopt
from itertools import groupby, count
from operator import itemgetter
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats

# init main
def main(argv):
    hairs_file = ''
    hapcut_file = ''
    bam_file = ''
    out_file = ''
    help = 'fake_hairs_from_table.py -r <chrom:start-end> -i <input.bam> -v <methyl.tsv> -o <output.hairs>'
    try:
        opts, args = getopt.getopt(argv,"r:i:v:o:",["region=", "input=", 'variant=', "output="])
    except getopt.GetoptError:
        print help
        sys.exit(2)
    for opt, arg in opts:
        if opt == '--help':
            print help
            sys.exit()
        elif opt in ("-r", "--region"):
		region = arg
		chrom = arg.split(':')[0]
		istart = int(arg.split(':')[1].split('-')[0])
		iend   = int(arg.split(':')[1].split('-')[1])
        elif opt in ("-i", "--input"):
            bam_file = arg
        elif opt in ("-o", "--output"):
            out_file = arg
	elif opt in ("-v", "--variant"):
            var_file = arg
        else:
            assert False, "unhandled option"
            
    out_fp = open(out_file, 'w')
    bam_fp = pysam.AlignmentFile(bam_file, "rb")
    mergedvcf = pd.read_csv(var_file, header=None, sep='\t')
	
    sys.stderr.write("Sucessful init: fake_hairs_from_table.py, starting...\n")
    sys.stderr.flush()
    fakeHairs(bam_fp, mergedvcf, out_fp, chrom, istart, iend)
    out_fp.close()
    bam_fp.close()

def hamming_dist(str1, str2):
    difference = 0
    for x,y in zip(str1, str2):
        if x != y:
            difference += 1
    return difference

def get_refmap(bamread):
    bpos = bamread.get_aligned_pairs(matches_only=True)
    positions = [i[0] for i in bpos]
    refpos =  [i[1] for i in bpos] # positions in reference
    refmap = dict(zip(refpos, positions))
    return refmap

def get_matched_bases_in_read(bamread, in_pos, refmap):
    outseq = [bamread.seq[refmap[i]] if i in refmap else 'N' for i in in_pos]
    outqual = [bamread.qual[refmap[i]] if i in refmap else '.' for i in in_pos]
    outseq = ''.join(outseq)
    outqual = ''.join(outqual)
    return outseq, outqual

def get_position_in_read(bamread, in_pos, refmap):
    outpos = None
    if in_pos in refmap:
        outpos = refmap[in_pos]
    return outpos

def reverse_compl(seq):
    translate = {'A':'T', 
                 'T':'A', 
                 'C':'G',
                 'G':'C', 
                 'N':'N'}
    return ''.join([translate[s] for s in seq])

def dummyScore(log_val, prior_mean_hap1, prior_mean_hap2):
	zscore1 = abs(log_val - prior_mean_hap1)/0.8 # dummy SD, fixed for all
	pval1 = stats.norm.sf(zscore1)
	zscore2 = abs(log_val - prior_mean_hap2)/0.8 # dummy SD, fixed for all
	pval2 = stats.norm.sf(zscore2)
	pval_hap1_nothap2 = (pval1 * (1-pval2))/(pval1*(1-pval2) + pval2*(1-pval1)) # p(hap1) * (1-p(hap2)) = P(hap1|not hap2)
	pval_hap2_nothap1 = (pval2 * (1-pval1))/(pval2*(1-pval1) + pval1*(1-pval2)) # p(hap1) * (1-p(hap2)) = P(hap1|not hap2)
	#print pval_hap1_nothap2, pval_hap2_nothap1
	if pval_hap2_nothap1 > pval_hap1_nothap2:
		return 1, toPhred(pval_hap1_nothap2)
	else:
		return 0, toPhred(pval_hap2_nothap1)
		
def toPhred(val):
    return chr(int(-10*np.log10(val)) + 33), int(-10*np.log10(val))

def fakeHairs(bam_fp, mergedvcf, out_new_hair, chrom, istart, iend):
	counter = 0
	for read in bam_fp.fetch(region=chrom, start=istart, end=iend):
		counter += 1
		if read.is_reverse:
			continue
		if counter % 100 == 0:
			out_new_hair.flush()
		refmap, ipd_values = None, None
		start, end = read.reference_start, read.reference_end
		read_variants = mergedvcf[(mergedvcf[1] >= start) & (mergedvcf[1] <= end)]
		if len(read_variants) != 0:
			refmap = get_refmap(read)
			ipd_values = [int(i) for i in read.get_tag('ip').strip('S').split(',')]
		phased_variants = []
		quals = ''
		for varid, varline in read_variants.iterrows():
			varid += 1 # this gives us the correct line numbers
			refpos = varline[1]
			allele = 0
			basepos = get_position_in_read(read, refpos, refmap)
			mean1 = varline[2]
			mean2 = varline[3]
			if mean1 > mean2: # swap so in order of size
				tmp = mean1
				mean1 = mean2
				mean2 = tmp
			if basepos==None: continue
			ipd_value = np.log(ipd_values[basepos])
			allele, qualval = dummyScore(ipd_value, mean1, mean2)
			if qualval[1] < 13 : # minimum quality value cutoff!
				continue
			phased_variants.append((varid, allele))
			# new quality score calculation code!
			quals += qualval[0]
		if len(phased_variants) < 2:
			continue # don't print nonuseful fragments!
		groups = groupby(phased_variants, key=lambda item, c=count():item[0]-next(c))
		tmp = [list(g) for k, g in groups]
		blockcount = len(tmp)
		outstr = ''
		outstr += str(blockcount) + " " + read.qname + " "
		for site in tmp:
			outstr += str(site[0][0]) + " "+''.join([str(y[1]) for y in site]) + " "
		outstr += quals
		print >>out_new_hair, outstr

# run the program if called from the command line
if __name__ == "__main__":
   main(sys.argv[1:])


