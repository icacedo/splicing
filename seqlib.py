# Assignment 1: Use argparse and libraries
# Never copy-paste functions from one program to another
# make a library instead and import functions
import gzip
import random
import math

def read_fasta(filename):
	
	name = None
	seqs = []
	
	fp = None
	if filename.endswith('.gz'): 
		fp = gzip.open(filename, 'rt')
	else: 
		fp = open(filename)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:	
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()
	
def rev_comp(seq):

	comp = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')
	anti = seq.translate(comp)[::-1]
	return anti
	
def read_fastq(filename):
	
	if filename.endswith('.gz'): 
		fp = gzip.open(filename, 'rt')
	else: 
		fp = open(filename)

	while True:
		name = fp.readline()
		seq = fp.readline()
		plus = fp.readline()
		qual = fp.readline()
		if name == '': break
		yield(name, seq, qual)
		
def random_dna(length, A, C, G, T):
	
	assert(math.isclose(1, A+C+G+T, abs_tol = 0.001))
	seq = []
	for i in range(length):
		r = random.random()
		if r < A:		seq.append('A')
		elif r < A + C:	seq.append('C')
		elif r < A + C +G:	seq.append('G')
		else:			seq.append('T')
	return ''.join(seq)
		
		

		
		
		
		
		
		
		
		
		
		
		
		
		
		

