import sys
import gzip
import itertools as it
import modelib as ml
import isoform as iso
import seqlib as sql

seqs = []
labels = []
for label, seq in sql.read_fasta(sys.argv[1]):
	seqs.append(seq)
	labels.append(label)

#### re-write of the api algorithm, with more parameters ######################
# avg exon size in C. elegans: 200.7bp
# median: 123bp, smallest: 7bp, largest 7569bp
# avg intron size: 47bp
# median: 65bp, only 3 <= 25bp
# largest introns: 100912bp and 21230bp
# use GU/AG or GT/AG rule for donor/acceptor sites
# are the sequences in data/sequences/ translated mRNA or cDNA?
# need to include min intron, min exon, max splice, flank

# test seq is 40bp, same 20 bp sequence twice
seqs = ['ATATATCGTCGATCAGCCGTATATATCGTCGATCAGCCGT']
def splice_sites(seqs):

	don_sites = []
	acc_sites = []
	seq_lens = []
	for i in range(len(seqs)):
		for j in range(len(seqs[i])-1):
			if seqs[i][j:j+2]=='GT':
				don_sites.append(j+1)
			if seqs[i][j:j+2]=='AG':
				acc_sites.append(j+2)
			else: continue
		seq_lens.append(len(seqs[i]))	
		
	return don_sites, acc_sites, seq_lens

# need to filter after/during/while using the API algorithm
# flanking regions are not a part of the exon

flank5 = 5
minintron = 8
minexon = 7
sites = splice_sites(seqs)
d_sites = sites[0]
a_sites = sites[1]
s_lens = sites[2]

#print(d_sites, a_sites)

for i in range(1,len(d_sites)):
	for d in it.combinations(d_sites,i):
		for a in it.combinations(a_sites,i):
			# d and a are of equal length
			# can use same index to get interior exons
			#print(d,a)
			for don, acc in zip(d,a):
				# remove splice site if the donor comes after the acceptor
				if don >= acc: continue
				# remove short introns
				i_len = acc-don
				if i_len < minintron: continue
				# remove short 5' exon
				e_len5 = don - flank5
				if e_len5 < minexon: continue
				# remove short 3' exon
				e_len3 = s_lens[i-2] - acc
				if e_len3 < minexon: continue
			for j in range(len(d)):	
				print(j)
				print(d[j], a[j])				
			
				
				'''
				if i_len < minintron: continue
				e_len = 
				'''
				# think
				# if there are only one donor and acceptor sites
				# the exons will be on the flanks
				# does the min exon size
				# work for only one of the exons on the flank?
				
			
		






















	
