import os, sys


class BlockVariant:
	def __init__ (self, variantline):
    	# variant_id haplotype_1 haplotype_2 chromosome position refallele variantallele genotype allele_counts:genotype_likelihoods:delta:MEC_variant
		ll = variantline.strip().split("\t")
        var_id, hap1, hap2, chrom, pos, r_allele, v_allele, genotype, info_str = ll
        self.var_id, self.hap1, self.hap2, self.pos = int(var_id), int(hap1), int(hap2), int(pos)
        allele_counts, genotype_likelihoods, delta, MEC_variant = info_str.split(":")
        self.ref_count, self.alt_count = map(int, allele_counts.split(","))
        gen_00, gen_01, gen_11 = map(float, genotype_likelihoods.split(","))
        self.gen_like = {"0/0":gen_00, "0/1":gen_01, "1/1":gen_11}
		self.delta = float(delta)
		self.MEC_variant = MEC_variant


class Block:
	def __init__ (self, blockline):
        # "BLOCK: offset:" first_variant_block "len:" length_of_block "phased": phased_variants_block SPAN: lengthspanned MECscore score fragments #fragments  
        ll               = blockline.strip().split()
		self.offset      = int(ll[2])
		self.total_len   = int(ll[4])
		self.phased      = int(ll[6])
		self.span        = int(ll[8])
		self.MECscore    = float(ll[10])
		self.fragments   = int(ll[12])
		variants = [] # default to empty

	def addVariant (variantline):
       	variants.append(BlockVariant(variantline))


class HapCutReader:
    
    def __init__ ( self, fn ):
        self.blocks = read_file_to_blocks (fn)
        

    def read_file_to_blocks (self, fn):
        with open (hapcutoutfn) as f:
            f.readline()
            blocks = []
			currBlock = None
            prevBlock = False
            for l in f:
                if l[0] == "B":
                	if prevBlock:
	                    yield block, snpDict
	                else:
	                	prevBlock = True
	                currBlock = Block(l)
                else:
                	currBlock.addVariant(l)

                        

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
		positions       = []                   # an array with the indices covered by read
		alleles         = []                   # a matched array with the allele calls at each position
		for i in range(2, len(hairlist)-1, 2): 
			positions.append(hairlist[i])
			positions.append(hairlist[i+1])
		self.qualities  = hairlist[-1]         # a matched arary of the qualities of allele calls

class HairReader:

	def __init__ (self, fn):
		self.reads = []
		with open (fn) as f:
			for l in f:
				self.reads.append(HapCutRead(l))


