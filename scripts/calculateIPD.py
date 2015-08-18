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

	helpline = 'calculateIPD.py -i <in.bam> -r <ref.fa> -c <chrom:start-end> -p <prefix>'
	try:
		opts, args = getopt.getopt(argv,"i:r:c:p:",["in=","ref=", "region=", "prefix="])
	except getopt.GetoptError:
		print helpline
		sys.exit(2)
	for opt, arg in opts:
		if opt == '--help':
			print helpline
			sys.exit()
		elif opt in ("-i", "--in"):
			bamfile = arg
		elif opt in ("-r", "--ref"):
			reffile = arg
		elif opt in ("-c", "--region"):
			region = arg
			chrom = arg.split(':')[0]
			start = int(arg.split(':')[1].split('-')[0])
			end   = int(arg.split(':')[1].split('-')[1])
		elif opt in ("-p", "--prefix"):
			positionfile = arg
			motiffile = arg
		else:
			assert False, "unhandled option"

	assert pysam.Samfile(bamfile, 'rb'), 'ERROR: Cannot open bam file for reading.'
	assert open(bamfile + '.bai', 'rb'), 'ERROR: bam file is not indexed!'
	bam_fp = pysam.Samfile(bamfile, 'rb')    

	motiffile = motiffile + "_priors_" + region + ".tsv"
	positionfile = positionfile + "_positions_" + region + ".tsv"

	assert open(motiffile, 'wb'), 'ERROR: Cannot open output file.'
	motif_fp = open(motiffile, 'wb')
	assert open(positionfile, 'wb'), 'ERROR: Cannot open output file.'
	pos_fp = open(positionfile, 'wb')

	try:
		sequences = dict( (s.name, s) for s in ht.FastaReader(reffile) )
	except:
		sys.stderr.write("ERROR: Cannot open reference file.\n")
		sys.exit(2)

	if chrom not in sequences:
		sys.stderr.write("ERROR: Region not in reference file.\n")
		sys.exit(2)

	sys.stderr.write("Sucessful init: calculateIPD.py, starting...\n")
	sys.stderr.write("chrom: %s start: %s end: %s\n" % (chrom, str(start), str(end)))
	sys.stderr.flush()
	getIPDfromBAM(bam_fp, sequences, chrom, start, end, motif_fp, pos_fp)
	motif_fp.close()
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
    
def getIPDfromBAM(bam_fp, sequences, chrom, start, end, motif_fp, pos_fp):
	ipd_by_position = dict()
	motif_ipd_values_read = dict()
	count=0
	for read in bam_fp.fetch(reference=chrom, start=start, end=end):
		count += 1
		if read.is_reverse:
			continue
		if count % 100 == 0:
			sys.stderr.write("Read %s reads from file\n" % str(count))
		ipd_values = [int(i) for i in read.get_tag('ip').strip('S').split(',')]
		rstart, rend = read.pos, read.aend
		rlen = rend-rstart
		refseq = str(sequences[chrom][rstart:rend])
		refmap = get_refmap(read)
		for refpos in range(rstart+6,rend-6):
			if (refpos < start) | (refpos > end): continue # don't duplicate effort
			basepos = get_position_in_read(read, refpos, refmap) # translate the refpos value to a read-based coord
			if basepos == None: continue
			motif = refseq[refpos-5-(rstart):refpos+6-rstart]
			ref_3mer = motif[4:7]
			read_3mer = get_matched_bases_in_read(read, range(refpos-1,refpos+2), refmap)
			if ref_3mer != read_3mer: continue
			ipd_value = ipd_values[basepos] 
			if refpos in ipd_by_position:
				ipd_by_position[refpos].append(ipd_value)
			else:
				ipd_by_position[refpos] = [ipd_value]
			if motif in motif_ipd_values_read:
				motif_ipd_values_read[motif].append(ipd_value)
			else:
				motif_ipd_values_read[motif] = [ipd_value]
	sys.stderr.write("Finished calculating IPDs. Writing positional output\n")
	pos_fp.write("#ref_chrom\tref_pos\tipd_at_position\n")
	for key in sorted(ipd_by_position.keys()):
		pos_fp.write('\t'.join([str(x) for x in [chrom, key, ','.join([str(a) for a in ipd_by_position[key]])]]) + '\n')
	pos_fp.flush()
	sys.stderr.write("Finished writing positional data. Writing priors to file\n")
	motif_fp.write("#11mer_motif\tipd_at_mid_pos\n")
	for key in motif_ipd_values_read.keys():
		motif_fp.write('\t'.join([key, ','.join([str(a) for a in motif_ipd_values_read[key]])]) + '\n')
	motif_fp.flush()
	sys.stderr.write("Successfully completed.\n")
	sys.stderr.flush()
	return 0

# run the program if called from the command line
if __name__ == "__main__":
   main(sys.argv[1:])

