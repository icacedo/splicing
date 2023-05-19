# code exon/intron markov model
import isoform
import sys
import math

exin_seqs = []
with open(sys.argv[1], 'r') as fp:	
	exin_seq = ''
	for line in fp.readlines():
		line = line.rstrip()
		exin_seq = line
		exin_seqs.append(exin_seq)


seqs = [
	'CTTGGACATATCGCATATCGATTTGCTATTCCGGATGCTAATGAGCTCCGCCATTACGCCGCTGGATCATCAACC', 
	'TTGCCAACCAACGGAGTTGTCGCAGACGAAAGAAAGTATACATCTACAAGCTGAACGACAGGCTCACACAGTGCGGAATG', 
	'ATGTCTCTTAAGATTTCTTATTTTCGAACACGTTGCACATCCGCAGAACCTTGAAAGAGCTTCTGGAAGT'
]
s = 'ATGGGCGCGTTATATGGCTACGATAATATCGACTATC'
#s = 'ACTGACTGAC'
order = 3


context = {}
for seq in exin_seqs:
	for i in range(len(seq)-order):
		prev = seq[i:i+order]
		now = seq[i+order]
		if prev not in context:
			context[prev] = now
		else:
			context[seq[i:i+order]] += now
#print(context)

# this is correct compared to exon.mm
for nts in sorted(context):
	#print(nts, context[nts])
	print(nts)
	A = 0
	C = 0
	G = 0
	T = 0
	for nt in context[nts]:
		if nt == 'A': A += 1
		if nt == 'C': C += 1
		if nt == 'G': G += 1
		if nt == 'T': T += 1
		d = int(len(context[nts]))
	print(A/d, C/d, G/d, T/d)





'''
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

print(math.log(allk['AAAA']/allk_counts, 2))

kmers, total = kcounts(s, k)

for i in kmers:
	print(i, kmers[i]/total)

print('********')
# idk how to work
#print(isoform.create_markov(sys.argv[1], 2, 0, 1111))
'''





