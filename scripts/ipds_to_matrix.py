#!/usr/bin/env python

import pysam
import HTSeq as ht
import math
import os, sys, getopt

# init main
def main(argv):
    hairs_file = ''
    hapcut_file = ''
    bam_file = ''
    out_file = ''
    help = 'ipds_to_matrix.py -r <chr:start-stop> -i <input.bam> -o <output.tsv> -f <ref.fa>'
    try:
        opts, args = getopt.getopt(argv,"r:i:o:f:",["region=", "input=", "output=", "ref="])
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
		elif opt in ("-f", "--ref"):
			ref_file = arg
		else:
			assert False, "unhandled option"
            
    out_fp = open(out_file, 'w')
    bam_fp = pysam.AlignmentFile(bam_file, "rb")
    sequences = dict( (s.name, s) for s in ht.FastaReader(ref_file) )
	
    sys.stderr.write("Sucessful init: ipds_to_matrix, starting...\n")
    sys.stderr.flush()
    ipds_to_matrix(bam_fp, out_fp, sequences, chrom, istart, iend)
    out_fp.close()
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

def ipds_to_matrix(bam_fp, out_fp, sequences, chrom, istart, iend):
	counter = 0
	for read in bam_fp.fetch(region=chrom, start=istart, end=iend):
		counter += 1
		if counter % 100 == 0:
			out_fp.flush()
		readname = read.qname
		refmap, ipd_values = None, None
		start, end = read.reference_start, read.reference_end
		refseq = str(sequences[chrom][start:end])
		refmap = get_refmap(read)
		ipd_values = [int(i) for i in read.get_tag('ip').strip('S').split(',')]
		for refpos in refmap:
			ref_3mer = refseq[refpos-start-1:refpos-start+2]
			read_3mer = get_matched_bases_in_read(read, range(refpos-1,refpos+2), refmap)
			if ref_3mer != read_3mer: continue
			out_fp.write('\t'.join([readname, str(refpos), str(ipd_values[refmap[refpos]])]) + '\n')
			
# run the program if called from the command line
if __name__ == "__main__":
   main(sys.argv[1:])


