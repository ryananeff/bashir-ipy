#!/usr/bin/env python

# imports
import os, sys, getopt
import pysam
from itertools import groupby
import pandas as pd
import numpy as np
global count
count = 0

# init main
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
    assert open(bam_file + '.bai', 'rb'), 'ERROR: bam file is not indexed!'
    bam_fp = pysam.Samfile(bam_file, 'rb')

    if out_file==None:
        out_file = bam_file + ".ann_haplotypes_" + time.strftime("%m%d%y_%H%M%S") + '.bam'

    assert pysam.AlignmentFile(out_file, "wb", template=bam_fp), 'ERROR: Cannot open output file for writing.'
    out_fp = pysam.AlignmentFile(out_file, "wb", template=bam_fp)

    assert open(hairs_file), 'ERROR: Cannot access hairs file.'
    assert open(hapcut_file), 'ERROR: Cannot open hapcut file.'

    hair_reader = HairReader(hairs_file)
    block_reader = HapCutReader(hapcut_file)
    stats_file = out_file + ".interblock_stats.tsv"
    sys.stdout.write("Loaded greedy_partitoner.py, beginning execution. \n")
    sys.stdout.flush()
    tag_reads(bam_fp, hair_reader, block_reader, out_fp) #begin tagging reads
    interblock_stats(hair_reader, block_reader, stats_file) #generate interblock stats
    bam_fp.close()
    out_fp.close()

# end of main

### CLASSES ###

class BlockVariant:
    def __init__ (self, variantline):
        # variant_id haplotype_1 haplotype_2 chromosome position refallele variantallele genotype allele_counts:genotype_likelihoods:delta:MEC_variant
        ll = variantline.strip().split("\t")
        var_id, hap1, hap2, chrom, pos, r_allele, v_allele, genotype, info_str = ll
        self.chrom, self.r_allele, self.v_allele, self.info_str = chrom, r_allele, v_allele, info_str
        self.var_id, self.hap1, self.hap2, self.pos = int(var_id), hap1, hap2, int(pos)
        allele_counts, genotype_likelihoods, delta, MEC_variant = info_str.split(":")[0:4]
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
        self.variants 	 = dict() # default to empty
        self.variant_ids = set()
        self.chrom = 0
        self.start = 0
        self.end = 0
        self.informative_reads = []
        self.read_count = 0
        self.read_set = frozenset()

    def __repr__ (self):
        return "<Block, offset_id: %s>" % str(self.offset)

    def addVariant(self, variantline):
        variant = BlockVariant(variantline)
        self.variants[variant.var_id] = variant
        self.variant_ids.add(variant.var_id)
        self.updatePosition()

    def updatePosition(self): # we need to do this because sometimes the variant isn't associated with a block
        positions = []
        chrom = None
        for k,variant in self.variants.iteritems():
            if chrom == None:
                chrom = variant.chrom
            positions.append(variant.pos)
        self.chrom = chrom
        self.start = np.min(positions)
        self.end = np.max(positions)

    def addReadsToBlock(self, read_dict):
        self.informative_reads = []
        for k,read in read_dict.iteritems():
            read_ids = [var[0] for block in read.blocks for var in block]
            if len(set(read_ids).intersection(set(self.variant_ids))) > 0:
                self.informative_reads.append(read)
        self.read_count = len(self.informative_reads)
        self.read_set = frozenset([x.read_id for x in self.informative_reads])

    def concordance(self, input_reads):
        ''' this should return a dict of (#T,#F) tuples per variant
         each element is a variant's concordance with the reads
         using the read's haplotype information, we can establish whether the read's phasing
         is consistent with how the variant was phased '''
        variant_concord = dict()
        support_reads_hap2 = 0
        against_reads_hap2 = 0
        support_reads_hap1 = 0
        against_reads_hap1 = 0
        for k,variant in self.variants.iteritems():
            for read in input_reads:
                if variant.var_id in read.positions:
                    read_allele = read.alleles[read.positions.index(variant.var_id)]
                    hapstate = read.haplotypes[self.offset]
                    if hapstate == 2:
                        if variant.hap2 != read_allele:
                            against_reads_hap2 += 1
                        else:
                            support_reads_hap2 += 1
                    else:
                        if variant.hap1 != read_allele:
                            against_reads_hap1 += 1
                        else:
                            support_reads_hap1 += 1
        variant_concord[self.offset] = {"hap1": (support_reads_hap1, against_reads_hap1), 
                                               "hap2": (support_reads_hap2, against_reads_hap2)}
        return variant_concord

    def variant(self, var_id):
        try:
            return self.variants[var_id]
        except:
            return None

    def interblock_reads(self, input_reads):
        out_reads = []
        for read in input_reads:
            if read.read_id not in self.read_set:
                if read.chrom == self.chrom:
                    if ((read.start < self.end) & (read.end > self.end)) | \
                        ((read.end > self.start) & (read.start < self.start)) | \
                        ((read.end <= self.end) & (read.start >= self.start)):
                        out_reads.append(read)
        return out_reads

class HapCutReader:

    def __init__ ( self, fn ):
        self.fn = fn
        self.blocks = dict()
        self.translate = dict()
        for block in self.read_file_to_blocks(fn):
            self.blocks[block.offset] = block
            for v in block.variant_ids:
                self.translate[v] = block.offset

    def loc(self, block_id):
        try:
            return self.blocks[self.translate[block_id]]
        except:
            return None

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
        #Column 1 is the number of consecutive set of SNPs covered by the fragment, NOT haplotype blocks.
        #Column 2 is the fragment id. 
        #Column 3 is the offset of the first block of SNPs covered by the fragment followed by the alleles at the SNPs in this block.
        #Column 5 is the offset of the second block of SNPs covered by the fragment followed by the alleles at the SNPs in this block.
        #...
        #The last column is a string with the quality values (Sanger fastq format) for all the alleles covered by the fragment (concatenated for all blocks). 
        #For example, if a read/fragment covers SNPs 2,3 and 5 with the alleles 0, 1 and 0 respectively, then the input will be:
        #2 read_id 2 01 5 0 AAC
        #Here AAC is the string corresponding to the quality values at the three alleles. The encoding of 0/1 is arbitrary but following the VCF format, 0 is reference and 1 is alternate. 
        hairlist = hairline.strip().split()
        self.blockcount = 0                # this information must be determined afterwards 
        self.read_id    = hairlist[1]      # read_id
        self.blocks     = []		       # an array of tuples corresponding to blocks
        self.positions  = []
        self.alleles    = []
        self.chrom      = None
        self.start      = None
        self.end        = None
        self.haplotypes = dict()             # an array of {"block_offset":"haplotype"} 
                                           # after partitioning
        for i in range(2, len(hairlist)-1, 2):
            position = int(hairlist[i])
            allele = hairlist[i+1]
            block = zip(range(position, position+len(allele)), allele)
            self.blocks.append(block)
            self.positions.extend(range(position, position+len(allele)))
            self.alleles.extend(allele)
            self.qualities  = hairlist[-1]         # a matched arary of the qualities of allele calls

    def __repr__(self):
        return "<HapCutRead, read_id: %s>" % str(self.read_id)

    def haplotype_fields(self):
        haps = ";".join([','.join([str(key), str(self.haplotypes[key])]) for key in self.haplotypes])
        haptag = [("ZH", haps), ("ZB", int(self.blockcount))]
        return haptag
    
    def addGenomicPositions(self, block_reader):
        arr = []
        chrom = None
        for position in self.positions:
            b = block_reader.loc(position)
            if b == None:
                continue
            if chrom == None:
                chrom = b.chrom
            arr.append(b.variant(position).pos)
        if len(arr) > 0:
            self.chrom = chrom
            self.start = np.min(arr)
            self.end = np.max(arr)
        else:
            self.chrom = '*'
            self.start = None
            self.end = None
    
class HairReader:

    def __init__ (self, fn):
        self.fn = fn
        self.reads = dict()
        with open (fn) as f:
            for l in f:
                newread = HapCutRead(l)
                self.reads[newread.read_id] = newread
        self.read_set = frozenset(self.reads.keys())

    def __repr__ (self):
        return "<HairReader, filename: %s>" % self.fn

    def loc(self, read_id):
        try:
            return self.reads[read_id]
        except:
            return None

### FUNCTIONS ###

'''

tag_reads()

Usage: Tags reads from a bam file corresponding to a particular haplotype, with haplotype
definitions from HapCut, under the optional tag "ZH".

Inputs:
    bam_fp
    A pysam.Samfile object pointing to the input file
    hair_reader
    A HairReader object pointing to the hairs file.
    block_reader
    A HapcutReader object pointing to the hapcut file.
    out_fp
    A pysam.AlignmentFile pointing to the output bam.
Outputs:
    (none - writes to out_fp)

'''

def tag_reads(bam_fp, hair_reader, block_reader, out_fp):
    ''' tag_reads(bam_fp, hair_reader, block_reader, out_fp)'''
    global count
    for bamread in bam_fp.fetch():
        count += 1
        if (count % 100) == 0:
            sys.stdout.write("\rWritten %s lines to output." % str(count))
            sys.stdout.flush()
        if bamread.query_name in hair_reader.read_set:
            read = hair_reader.loc(bamread.query_name)
            read = greedy_partition(read, block_reader)
            bamread.tags += read.haplotype_fields()      # add the haplotype information
        out_fp.write(bamread)

'''
greedy_partition()
Ryan Neff

inputs:
read
    a HapCutRead object
block_reader
    of the type HapCutReader

outputs:
    the original read, now with haplotype information.

translate hairfile alleles into blockvar IDs
get alleles in each read spanning blockvars
determine alleles for the two blocks from blockvar
partition read based on locally most probable alignment

'''

def greedy_partition(read, block_reader):
    # it turns out that the blocks provided in a hapcut file don't actually correspond to real blocks
    # just contiguous alleles?
    read.blocks = [] # because they are useless
    lastBlock = -1
    for ix, pos in enumerate(read.positions):
        currBlock = block_reader.loc(pos)
        if currBlock == None:
            continue
        currBlock = currBlock.offset
        if currBlock != lastBlock:
            read.blocks.append([])
        read.blocks[-1].append((pos, read.alleles[ix]))
        lastBlock = currBlock
    read.blockcount = len(read.blocks)
    for readblock in read.blocks:
        positions = [x[0] for x in readblock]
        alleles = [x[1] for x in readblock]
        allele_state = []
        offset = positions[0]
        hap = 0
        block = block_reader.loc(offset) #retrieve block the read is in
        offset = block.offset
        if block == None:
            continue
        for ix, varpos in enumerate(positions):
            blockvar = block.variant(varpos)
            if blockvar.hap1 == alleles[ix]:
                allele_state.append(-1)
            elif blockvar.hap2 == alleles[ix]:
                allele_state.append(1)
            else:
                #sys.stderr.write("\nWarning: read allele matched no haplotypes.")
                #sys.stderr.write("\nHair read: %s" % read.read_id)
                #sys.stderr.write("\nAlleles: %s" % str(alleles[ix]))
                #sys.stderr.write("\n Hap 1: %s, Hap 2: %s\n" % (blockvar.hap1, blockvar.hap2))
                sys.stderr.flush()
                continue
        if len(allele_state) < 1:
            sys.stderr.write("Warning: no haplotype information in read.\n")
            sys.stderr.write("\nHair read: %s" % read.read_id)
            sys.stderr.flush()
            hap = -1
        if sum(allele_state) < 0:
            hap = 1
        elif sum(allele_state) > 0:
            hap = 2
        else:
            hap = 0
        read.haplotypes[offset] = hap
    return read

'''
interblock_stats()

Usage: Creates a tab-separated values file with statistics about reads overlapping
between nearby blocks, and finds the concordance of these interblock reads
with haplotypes in other blocks. 

inputs:
    hair_reader
        A HairReader object
    block_reader
        A HapcutReader object
    out_stats
        A string where the .tsv should be written. Defaults
        to the hairs filename given in the input + 'interblock_stats.tsv'
outputs:
    none-writes to file directly

'''

def interblock_stats(hair_reader, block_reader, out_stats):
    blockdist = []
    lastChr = None
    lastPos = None
    lastBlock = None
    lastReads = set()
    for k,read in hair_reader.reads.iteritems():
        if read.haplotypes == dict():
            read = greedy_partition(read, block_reader)
        read.addGenomicPositions(block_reader)
    for ix, key in enumerate(sorted(block_reader.blocks.keys())):
        sys.stdout.write('\r%s percent done.' % round(ix/float(len(block_reader.blocks))*100))
        sys.stdout.flush()
        block = block_reader.blocks[key]
        if block.read_set == set():
            block.addReadsToBlock(hair_reader.reads)
        currBlock = block.offset
        currChr = block.chrom
        currPos = block.start
        if lastBlock != None:
            if lastChr == currChr:
                interblock_reads = block.interblock_reads(lastBlock_obj.informative_reads)
                row=[lastBlock, currBlock, currChr, lastPos, currPos, currPos-lastPos,
                     lastBlock_obj.end-lastBlock_obj.start, block.end-block.start,
                     len(lastBlock_obj.variant_ids), len(block.variant_ids), 
                     len(lastBlock_obj.informative_reads), len(block.informative_reads),
                     len(list(bam_fp.fetch(region=lastChr + ':' + str(lastBlock_obj.start) + '-' + str(lastBlock_obj.end)))),
                     len(list(bam_fp.fetch(region=lastChr + ':' + str(block.start) + '-' + str(block.end)))),
                     len(interblock_reads),
                     len(list(bam_fp.fetch(region=lastChr + ':' + str(lastBlock_obj.end) + '-' + str(block.start)))),
                     lastBlock_obj.concordance(lastBlock_obj.informative_reads), 
                     block.concordance(block.informative_reads)]
                blockdist.append(row)
            else:
                continue
        lastBlock = currBlock
        lastBlock_obj = block
        lastChr = currChr
        lastPos = block.end
    header = ['block1', 'block2', 'chrom', 'block1_end', 'block2_start', 
              'distance', 'block1_size', 'block2_size', 'block1_variants', 'block2_variants', 
              'block1_informative_reads', 'block2_informative_reads', 'block1_reads', 'block2_reads',
              'informative_interblock_reads', 'all_interblock_reads', 'block1_concordance', 'block2_concordance']
    info = pd.DataFrame(blockdist, columns=header)
    info.to_csv(out_stats, sep="\t")

# run the program if called from the command line
if __name__ == "__main__":
   main(sys.argv[1:])

