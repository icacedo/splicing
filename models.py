import sys
import gzip
import itertools as it
import modelib as ml
import isoform as iso
import seqlib as sql
'''
seqs = []
labels = []
for label, seq in sql.read_fasta(sys.argv[1]):
	seqs.append(seq)
	labels.append(label)
'''
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
#########1######8######15##19#######28#####35##39#####
#########1######GT#####AG##GT#######GT#####AG##GT#####
seqs = ['ATATATCGTCGATCAGCCGTATATATCGTCGATCAGCCGT']
def splice_sites(seqs):

	seq_don_sites = []
	seq_acc_sites = []
	seq_lens = []
	for i in range(len(seqs)):
		d_sites = []
		a_sites = []
		for j in range(len(seqs[i])):
			if seqs[i][j:j+2]=='GT':
				d_sites.append(j+1)
			if seqs[i][j:j+2]=='AG':
				a_sites.append(j+1)
			else: continue
		seq_don_sites.append(d_sites)
		seq_acc_sites.append(a_sites)
		seq_lens.append(len(seqs[i]))	
		
	return seq_don_sites, seq_acc_sites, seq_lens

# need to filter after/during/while using the API algorithm
# flanking regions are not a part of the exon

# should flanking region be the same as minimum exon length?
flank5 = 5
flank3 = 5
minintron = 4
minexon = 4

sites = splice_sites(seqs)
don_sites = sites[0]
acc_sites = sites[1]
s_lens = sites[2]

isoforms = []
for i in range(len(s_lens)):
	for l in range(1,len(don_sites[i])+1):
		for d in it.combinations(don_sites[i],l):
			for a in it.combinations(acc_sites[i],l):
			
				# d and a are of equal length
				# why are unequal values of l ignored?
				# it.combinations() returns nothing if r > len(iterator)
				# if a don/acc list larger than the other, those iterations
				# are probably still done in the background
				# will this affect performance? idk yet
				# can use same index to get interior exons
				
				# remove short 5' exon
				if d[0] < flank5: continue
				
				# remove short 3' exon
				e_len3 = s_lens[i] - a[-1]
				if e_len3 < flank3: continue
				
				isof = []
				for don, acc in zip(d,a):
					# remove splice site if the donor comes after the acceptor
					if don >= acc: continue
					# remove short (includes interior) introns
					i_len = acc-don
					if i_len < minintron: continue
					isof.append((don,acc))
				# check for repeats
				if isof not in isoforms and isof != []:
					isoforms.append(isof)
print(isoforms)
		
			
	
print('****************')
seq='ATATATCGTCGATCAGCCGTATATATCGTCGATCAGCCGT'
txs, info = iso.all_possible(seq, 1, 1, 100, 1, gff=None)
for i in txs:
	print(i)



















	
