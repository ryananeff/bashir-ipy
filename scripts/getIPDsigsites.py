#!/usr/bin/env python

import pysam
import HTSeq as ht
import math
import os, sys, getopt

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

	helpline = 'getIPDsigsites.py -i <ipd_positions> -o <out_sig> -p <priors> -r <ref.fa>'
	try:
		opts, args = getopt.getopt(argv,"i:o:p:r:",["in=","out=", "prior=", "ref="])
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
		elif opt in ("-p", "--prior"):
			prior_file = arg
		elif opt in ("-r", "--ref"):
			ref_file = arg
		else:
			assert False, "unhandled option"

	pos_in = open(pos_in_file, 'r')
	out_sig = open(out_file, 'wb')
	prior_fp = open(prior_file, 'r')
	sequences = dict( (s.name, s) for s in ht.FastaReader(reffile) )
	
	sys.stderr.write("Sucessful init: getIPDsigsites.py, starting...\n")
	sys.stderr.flush()
	sys.stderr.write("
	getIPDfromBAM(bam_fp, sequences, chrom, start, end, pos_fp)
	pos_fp.close()
	bam_fp.close()


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
print >>sys.stderr, "init"
pos_in = open(sig_positions, 'rb') #'/hpc/users/neffr01/jason_new/hapcut_outputs/hg002_re_000000F/hapcut_qv13_mq10/hg002_qv13_redo_merged_haps.ipd.tsv'
block_reader = HapCutReader("/hpc/users/neffr01/jason_new/hapcut_outputs/hg002_re_000000F/hapcut_qv13_mq10/hg002_hapcut_000000F.hapcut")
bam_fp = pysam.AlignmentFile("/hpc/users/neffr01/jason_new/hapcut_outputs/hg002_re_000000F/hg002_000000F.new.merged.bam.rg.bam")
sequences = dict( (s.name, s) for s in ht.FastaReader("/hpc/users/neffr01/jason_new/contig_000000F.fa") )
#motif dict load?

print >>sys.stderr, "1. get set of significant IPD positions"
siteslist = [] # list of methylated positions
for line in pos_in:
    if (line[0] == "#") | (line[0] == "="):
        continue
    chrom, pos, ipdl_hap1, ipdl_hap2 = line.strip('\n').split('\t')
    pos = int(pos)
    curr_motif = str(sequences['000000F'][pos-5:pos+6])
    if curr_motif not in motif_dict:
        continue
    if (ipdl_hap1 == '')|(ipdl_hap2 == ''):
        continue
    ipd_hap1, ipd_hap2 = [int(x) for x in ipdl_hap1.split(',')], [int(x) for x in ipdl_hap2.split(',')]
    if min([ipd_hap1,ipd_hap2]) < 8:
        continue
    prior_ipds = [int(x) for x in motif_dict[curr_motif].split(',', 2001)[0:2000]]
    kstest_haps = stats.ks_2samp(ipd_hap1, ipd_hap2)
    kstest_prior1, kstest_prior2 = stats.ks_2samp(prior_ipds, ipd_hap1), stats.ks_2samp(prior_ipds, ipd_hap2)
    if (kstest_haps[1] > 0.01 ) | (min([kstest_prior1[1], kstest_prior2[1]]) > 0.01):
        continue
    logx, halflogx = np.linspace(0,10,100), np.linspace(0,10,50)
    kde_prior = kde_scipy(np.log(prior_ipds), logx)
    kde_hap1, kde_hap2 = kde_scipy(np.log(ipd_hap1), logx, bandwidth=0.25), kde_scipy(np.log(ipd_hap2), logx, bandwidth=0.25)
    obs1, obs2 = np.array([[np.log(i)] for i in ipd_hap1]), np.array([[np.log(i)] for i in ipd_hap2])
    obs = np.concatenate((obs1, obs2))
    g, h = mixture.GMM(n_components=1), mixture.GMM(n_components=1)
    x, y = g.fit(obs1), h.fit(obs2)
    #print "000000F", pos, x.means_[0][0], y.means_[0][0]
    siteslist.append(["000000F", pos, x.means_[0][0], y.means_[0][0]])

# run the program if called from the command line
if __name__ == "__main__":
   main(sys.argv[1:])

