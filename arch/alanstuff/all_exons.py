import itertools
import random
import sys
from Bio import SeqIO


for seq_record in SeqIO.parse("../data/sequences/777.fa","fasta"):
	print(seq_record.id)
	print(repr(seq_record.seq))
	print(len(seq_record))
	print(seq_record.seq)
seq = seq_record.seq


donor = []
acceptor = []

for i in range (len(seq)-1):
	if (seq[i] == 'G' and seq[i+1] == 'T'):
		donor.append(i)
	if (seq[i] == 'A' and seq[i+1] == 'G'):
		acceptor.append(i)

print("donor sites: ", donor)
print("acceptor sites: ", acceptor)


