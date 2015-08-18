#!/usr/bin/env python

import sklearn
import pysam
import math
from IPython import display
from sklearn.neighbors import KernelDensity
from sklearn import mixture
from scipy.stats import gaussian_kde
from statsmodels.nonparametric.kde import KDEUnivariate
from statsmodels.nonparametric.kernel_density import KDEMultivariate
import itertools
from scipy import stats
import numpy as np
import os, sys, getopt


def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)

# init main
def main(argv):
	bamfile = ''
	reffile = ''
	region = ''
	chrom = ''
	start = ''
	end = ''
	motiffile = ''
	positionfile = ''
	sequences = dict()

	helpline = 'getIPDsigsites.py -i <ipd_positions> -o <out_sig>'
	try:
		opts, args = getopt.getopt(argv,"i:o:",["in=","out="])
	except getopt.GetoptError:
		print helpline
		sys.exit(2)
	for opt, arg in opts:
		if opt == '--help':
			print helpline
			sys.exit()
		elif opt in ("-i", "--in"):
			pos_in_file = arg
		elif opt in ("-o", "--out"):
			out_file = arg
		else:
			assert False, "unhandled option"

	pos_in = open(pos_in_file, 'r')
	out_sig = open(out_file, 'wb')
	
	sys.stderr.write("Sucessful init: getIPDsigsites.py, starting...\n")
	sys.stderr.flush()
	getIPDfromBAM(pos_in, out_sig)
	pos_in.close()
	out_sig.close()


def get_refmap(bamread):
    bpos = bamread.get_aligned_pairs(matches_only=True)
    positions = [i[0] for i in bpos]
    refpos =  [i[1] for i in bpos] # positions in reference
    refmap = dict(zip(refpos, positions))
    return refmap

def get_matched_bases_in_read(bamread, in_pos, refmap):
    outseq = [bamread.seq[refmap[i]] if i in refmap else 'N' for i in in_pos]
    outseq = ''.join(outseq)
    return outseq

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
    
'''
haplotype joining by methlyation


pseudocode:
1. get set of significant IPD positions along the contig (by partitioned haplotype)
2. get set of blocks from hapcut file
3. for each pair of blocks: 
    a. look at reads that overlap a few ipd sites in each block, "call" read based on which model it falls into
        i. how to determine confidence interval? probably bootstrapping
    b. see if there is a pattern between haplotypes in each block and called haplotypes
    c. continue for whole contig/genome
'''
def getIPDfromBAM(pos_in, out_sig):
	print >>sys.stderr, "init"
	print >>sys.stderr, "1. get set of significant IPD positions"
	for line in pos_in:
		if (line[0] == "#") | (line[0] == "="):
			continue
		chrom, pos, ipd_at_position = line.strip('\n').split('\t')
		pos = int(pos)
		if (ipd_at_position == ''):
			continue
		ipds = [int(x) for x in ipd_at_position.split(',')]
		if len(ipds) < 20:
			continue
		ipds = np.log(ipds)
		logx, halflogx = np.linspace(0,10,100), np.linspace(0,10,50)
		kde_curr = kde_scipy(ipds, logx, bandwidth=0.25)
		obs = np.array([[i] for i in ipds])
		g = mixture.GMM(n_components=2)
		x = g.fit(obs)
		two_aic, two_bic = x.aic(obs), x.bic(obs)
		two_means = np.reshape(x.means_, (1,2))[0]
		two_weights = x.weights_
		g = mixture.GMM(n_components=1)
		x = g.fit(obs)
		one_aic, one_bic = x.aic(obs), x.bic(obs)
		one_mean = x.means_
		one_weight = x.weights_
		if two_aic+10 > one_aic: # cutoff
			continue
		if max(two_weights) > 0.75: # cutoff
			continue
		outline = ["000000F", pos, two_means[0], two_means[1], two_weights[0], two_weights[1], two_aic, two_bic, one_aic, one_bic]
		out_sig.write('\t'.join([str(a) for a in outline]) + '\n')
		out_sig.flush() # write for every site

# run the program if called from the command line
if __name__ == "__main__":
   main(sys.argv[1:])

