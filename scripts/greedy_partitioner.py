#!/usr/bin/env python

# imports
import os, sys, getopt
import pysam
from itertools import groupby

# init
def main(argv):
    hairs_file = ''
    hapcut_file = ''
    bam_file = ''
    out_file = ''
    help = 'greedy_partitioner.py -h <input.hairs> -c <input.hapcut> -i <input.bam> -o <output.ann.bam>'
    try:
        opts, args = getopt.getopt(argv,"h:c:i:o:",["hairs=","hapcut=", "input=", "output="])
    except getopt.GetoptError:
        print help
        sys.exit(2)
    for opt, arg in opts:
        if opt == '--help':
            print help
            sys.exit()
        elif opt in ("-h", "--hairs"):
            hairs_file = arg
        elif opt in ("-c", "--hapcut"):
            hapcut_file = arg
        elif opt in ("-i", "--input"):
            bam_file = arg
        elif opt in ("-o", "--output"):
            out_file = arg
        else:
            assert False, "unhandled option"

    assert pysam.Samfile(bam_file, 'rb'), 'ERROR: Cannot open bam file for reading.'
    bam_fp = pysam.Samfile(bam_file, 'rb')

    if out_file==None:
        out_file = bam_file + ".ann_haplotypes_" + time.strftime("%m%d%y_%H%M%S") + '.bam'

    assert pysam.AlignmentFile(out_file, "wb", template=bam_fp), 'ERROR: Cannot open output file for writing.'
    out_fp = pysam.AlignmentFile(out_file, "wb", template=bam_fp)

    assert open(hairs_file), 'ERROR: Cannot access hairs file.'
    assert open(hapcut_file), 'ERROR: Cannot open hapcut file.'

    tag_reads(bam_fp, hairs_file, hapcut_file, out_fp) #begin tagging reads
    bam_fp.close()
    out_fp.close()
# end of main()

### CLASSES ###

class BlockVariant:
    def __init__ (self, variantline):
        # variant_id haplotype_1 haplotype_2 chromosome position refallele variantallele genotype allele_counts:genotype_likelihoods:delta:MEC_variant
        ll = variantline.strip().split("\t")
        var_id, hap1, hap2, chrom, pos, r_allele, v_allele, genotype, info_str = ll
        self.chrom, self.r_allele, self.v_allele, self.info_str = chrom, r_allele, v_allele, info_str
        self.var_id, self.hap1, self.hap2, self.pos = int(var_id), hap1, hap2, int(pos)
        allele_counts, genotype_likelihoods, delta, MEC_variant = info_str.split(":")
        self.ref_count, self.alt_count = map(int, allele_counts.split(","))
        gen_00, gen_01, gen_11 = map(float, genotype_likelihoods.split(","))
        self.gen_like = {"0/0":gen_00, "0/1":gen_01, "1/1":gen_11}
        self.delta = float(delta)
        self.MEC_variant = MEC_variant
    def __repr__ (self):
        return "<BlockVariant, var_id: %s>" % str(self.var_id)


class Block:
    def __init__ (self, blockline):
        # "BLOCK: offset:" first_variant_block "len:" length_of_block "phased": phased_variants_block SPAN: 
        # lengthspanned MECscore score fragments #fragments

        ll               = blockline.strip().split()
        self.offset      = int(ll[2])
        self.total_len   = int(ll[4])
        self.phased      = int(ll[6])
        self.span        = int(ll[8])
        self.MECscore    = float(ll[10])
        self.fragments   = int(ll[12])
        self.variants 	 = [] # default to empty
        self.variant_ids = []

    def __repr__ (self):
        return "<Block, offset_id: %s>" % str(self.offset)

    def addVariant(self, variantline):
        variant = BlockVariant(variantline)
        self.variants.append(variant)
        self.variant_ids.append(variant.var_id)


class HapCutReader:

    def __init__ ( self, fn ):
        self.fn = fn
        self.blocks = list(self.read_file_to_blocks(fn))


    def read_file_to_blocks(self, fn):
        with open(fn) as f:
            currBlock = None
            for l in f:
                if l[0] == "B": # starting a new block
                    currBlock = Block(l)
                elif l[0] == "*": # ending a block
                    yield currBlock
                else:
                    currBlock.addVariant(l)

    def __repr__ (self):
        return "<HapCutReader, filename: %s>" % self.fn                    

class HapCutRead:

    def __init__ (self, hairline):
        #Column 1 is the number of blocks (consecutive set of SNPs covered by the fragment). 
        #Column 2 is the fragment id. 
        #Column 3 is the offset of the first block of SNPs covered by the fragment followed by the alleles at the SNPs in this block.
        #Column 5 is the offset of the second block of SNPs covered by the fragment followed by the alleles at the SNPs in this block.
        #...
        #The last column is a string with the quality values (Sanger fastq format) for all the alleles covered by the fragment (concatenated for all blocks). 
        #For example, if a read/fragment covers SNPs 2,3 and 5 with the alleles 0, 1 and 0 respectively, then the input will be:
        #2 read_id 2 01 5 0 AAC
        #Here AAC is the string corresponding to the quality values at the three alleles. The encoding of 0/1 is arbitrary but following the VCF format, 0 is reference and 1 is alternate. 
        hairlist = hairline.strip().split()
        self.blockcount = int(hairlist[0])     # number of blocks
        self.read_id    = hairlist[1]          # read_id
        self.blocks     = []		       # an array of tuples corresponding to blocks
        self.haplotypes = []		       # an array of {"block_offset":"haplotype"} 
                                # after partitioning
        for i in range(2, len(hairlist)-1, 2):
            position = int(hairlist[i])
            allele = hairlist[i+1]
            block = zip(range(position, position+len(allele)), allele)
            self.blocks.append(block)
            self.qualities  = hairlist[-1]         # a matched arary of the qualities of allele calls

    def __repr__(self):
        return "<HapCutRead, read_id: %s>" % str(self.read_id)

class HairReader:

    def __init__ (self, fn):
        self.fn = fn
        self.reads = []
        with open (fn) as f:
            for l in f:
                self.reads.append(HapCutRead(l))

    def __repr__ (self):
        return "<HairReader, filename: %s>" % self.fn

### FUNCTIONS ###

'''

tag_reads
7/5/2015

Usage: Tags reads from a bam file corresponding to a particular haplotype, with haplotype
definitions from HapCut, under the optional tag "ZH".

Inputs:
    bam_fp (pysam.Samfile)
        The pysam.Samfile object corresponding to a bam file on which the haplotype cuts were generated.
    hairs_file (string)
        The filename of a hairs file from Hapcut.
    hapcut_file (string)
        The filename of a hapcut blocks file.
    out_fp (pysam.AlignmentFile)
        Should be a writable file pysam file pointer.
Outputs:
    (no returns - it writes to out_fp)

'''

def tag_reads(bam_fp, hairs_file, hapcut_file, out_fp):
    count = 0
    sys.stdout.write('Started reading from BAM')
    sys.stdout.flush()
    read_array = HairReader(hairs_file).reads
    read_set = frozenset([x.read_id for x in read_array])
    block_array = HapCutReader(hapcut_file).blocks
    # let's group the reads by genomic position / block offset to speed up reading from the BAM file

    for bamread in bam_fp.fetch():
	count += 1
	if (count % 10000) == 0:
		sys.stdout.write('\rWritten %s reads to BAM' % str(count))
		sys.stdout.flush()
        if bamread.query_name in read_set:
            read = next((x for x in read_array if bamread.query_name == x.read_id))
            for readblock in read.blocks:
                block_id = readblock[0][0]
                block = next((x for x in block_array if block_id in x.variant_ids), None) #retrieve block from block_array
                if block == None: # we don't have the block information for the read
                    continue
                blockvar = next((x for x in block.variants if block_id == x.var_id), None) # retrieve variant from block
                read = greedy_partition(read, block_array)
                haps = ";".join([','.join([str(val) for val in block]) for block in read.haplotypes])
                haptag = [("ZH", haps), ("ZB", int(read.blockcount))]
                bamread.tags += haptag      # add the haplotype information
        out_fp.write(bamread)
    sys.stdout.write('\nSucessfully completed.\n')
    sys.stdout.flush()

'''
greedy_partition()
Ryan Neff

Usage: Runs a greedy algorithm to determine which haplotype block(s) a read
belongs to.

inputs:
    read (HapCutRead)
        A read from the hairs file.
    block_array (list of Block)
        All of the blocks in the hapcut file.

outputs:
    out_read (HapCutRead)
        The original read with haplotype information
'''

def greedy_partition(read, block_array):

    out_read = read
    out_read.haplotypes = [] # for safety
    for readblock in out_read.blocks:
        positions = [x[0] for x in readblock]
        alleles = [x[1] for x in readblock]
        allele_state = []
        offset = positions[0]
        hap = 0
        block = next((x for x in block_array if offset in x.variant_ids), None) #retrieve block the read is in
        for ix, varpos in enumerate(positions):
            blockvar = next((x for x in block.variants if x.var_id == varpos), None)
            if blockvar.hap1 == alleles[ix]:
                allele_state.append(-1)
            elif blockvar.hap2 == alleles[ix]:
                allele_state.append(1)
            else:
                sys.stderr.write("ERROR: read allele matched no haplotypes.")
                raise
        if len(allele_state) < 1:
            sys.stderr.write("Warning: no haplotype information in read.")
            hap = -1
        if sum(allele_state) < 0:
            hap = 1
        elif sum(allele_state) > 0:
            hap = 2
        else:
            hap = 0
        out_read.haplotypes.append((offset, hap))
    return out_read

# run the program if called from the command line
if __name__ == "__main__":
   main(sys.argv[1:])
