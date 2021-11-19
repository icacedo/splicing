import sys
import gzip
import itertools
import modelib as ml
import numpy as np
import isoform as iso
import seqlib as sql

seqs = []
labels = []
for label, seq in sql.read_fasta(sys.argv[1]):
	seqs.append(seq)
	labels.append(label)

#### re-write of the api algorithm, with more parameters ######################
# avg exon size in C. elegans: 200.7bp
# median: 123bp, smalest: 7bp, largest 7569bp
# avg intron size: 47bp
# median: 65bp, only 3 <= 25bp
# largest introns: 100912bp and 21230bp
# use GU/AG or GT/AG rule for donor/acceptor sites
# are the sequences in data/sequences/ translated mRNA or cDNA?
# need to include min intron, min exon, max splice, flank

# test seq is 43bp
seqs = ['ATATATCGTCGATCAGTCGATCGGTACACCTGGAGTTCACATT']
def splice_sites(seqs):

	don_sites = []
	acc_sites = []
	for i in range(len(seqs)):
		for j in range(len(seqs[i])-1):
			if seqs[i][j:j+2]=='GT':
				don_sites.append(j+1)
			if seqs[i][j:j+2]=='AG':
				acc_sites.append(j+2)
			else: continue

	return don_sites, acc_sites

# need to filter after/during/while using the API algorithm

'''
for d in d_sites(seqs):
	for a in a_sites(seqs):
		print(d-a)
'''


#for i in itertools.combinations([1,2,3,4],2):
#	print(i)
























	
