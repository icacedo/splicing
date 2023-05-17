# code exon/intron markov model
import isoform
import sys


s = 'ATGGGCGCGTTATATGGCTACGATAATATCGACTATC'
#s = 'ACTGACTGAC'
n = 4

def kcounts(seq, n):

	kmers = {}
	total = 0
	for i in range(len(seq)-n+1):
		kmer = seq[i:i+n]  
		if kmer not in kmers:
			kmers[kmer] = 1
			total += 1
		else:
			kmers[kmer] += 1
			total += 1
	return kmers, total

print(kcounts(s, n))

kmers, total = kcounts(s, n)

exin_seqs = []
with open(sys.argv[1], 'r') as fp:	
	exin_seq = ''
	for line in fp.readlines():
		line = line.rstrip()
		exin_seq = line
		exin_seqs.append(exin_seq)

#print(exin_seqs)

allk = {}
allk_counts = 0
for seq in exin_seqs:
	kmers, total = kcounts(seq, n)
	for kmer in kmers:
		if kmer not in allk:
			allk[kmer] = 1
		else:
			allk[kmer] += 1
	allk_counts += total

print(allk)
print(allk_counts)
print(allk['AAAA']/allk_counts)
'''
kmers, total = kcounts(s, k)

for i in kmers:
	print(i, kmers[i]/total)
'''
print('********')
# idk how to work
#print(isoform.create_markov(sys.argv[1], 2, 0, 1111))






