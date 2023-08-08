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
	total_iseqs = 0
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split()
		if line[2] == 'intron':
			print(line)
			beg = int(line[3]) - 1 
			end = int(line[4])
			print(beg, end)
			iseq = seq[beg:end]
			print(iseq, len(iseq))
			total_iseqs += len(iseq)
	print(total_iseqs, len(seq))

def total_iseqs(gff, fasta):

	seq = None
	for seqid, seq in ml.read_fastas(fasta):
		seq = seq

	with open(gff, 'r') as fp:
		total_iseqs = 0
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split()
			if line[2] == 'intron':
				beg = int(line[3]) - 1 
				end = int(line[4])
				iseq = seq[beg:end]
				total_iseqs += len(iseq)
		return len(seq), total_iseqs

seq_len, iseqs_len = total_iseqs(gff, fasta)

print(seq_len, iseqs_len)











