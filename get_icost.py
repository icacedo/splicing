# log of intron sequences/ total sequences
# use gff to get intron seqs
# calculate one intron cost for every isoform

import sys
import modelib as ml

seq_alt = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'

gff = sys.argv[1]
fasta = sys.argv[2]

seq = None
for seqid, seq in ml.read_fastas(fasta):
	seq = seq
print(seq, len(seq))

with open(gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split()
		if line[2] == 'intron':
			print(line)
			beg = int(line[3]) - 1 
			end = int(line[4])
			print(beg, end)
			print(seq[beg:end])

