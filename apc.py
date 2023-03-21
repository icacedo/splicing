'''
>ch.9940 IV:4806597-4807419 Gene:WBGene00044483
CTTTAACTGTTTCTCTATTCGTTGTATATCGATTTTTATTTTGAAAATTATTCACTCTATCAAGAAGAAGCCGAAATTCC
GCCCACTGCATCTTTTTTCAGTTGTGCTTCTCCAAAAACCATTTTCAGATGTCACAGGATCCCACGATCTCCAATTTTTT
AAATGGAAAACCCATCATCGCACCTGTCGGGCCATCTTCTTTATTAACCAGGTTACTTTGAAGCAATTATTTAATTTTTA
AGAAAAATTAATTCGAAATTTAAGATTACAATCGTTCCTGCCTCAAATGGAAGAAGCCAATCGAAATATGCCGGAAAACG
CGTCAGCTGACGCTTTTCATATTGAAAATGTCGAGGATGACAGCTCGGACGACAGCGATTCGAATTCAGAATCGATTTCG
GAAGAGGTACCTGTCACAACTGACGAAAAGAGCACAGCGAGTAGTACACAGAGGATTGAAATTGATTTAGACGTGTTCAA
GGAGATGAACAGCGCAGTTGATGACAGAGATGTCGGAATTCAGAATGTCGAAGCTCTTCCAGAAGCATTTCAGTCAAAAG
GAGAGGATTCGAAACCGGAAATTACAAAGAAACCGTTGATTGAAGAATTATGAAATTTATGTTTTTTCTATTGTTGTGAA
TTTTTTGTTGAATTTGTTACTTTTTGACCTACCCCTCGAAACATTTTCTCACATTGCGAAATGAATCAATAAAATTTCGA
TTTTAGATTTTTATTCATTTATCTCATTGGTTTCTATCCATTTTTGAATTACCTTCCCATTTTTATTTATCTCTGCCTTT
TCTAAAATTCTCCGAAAACTTTT
'''

seq3 = 'CTTTAACTGTTTCTCTATTCAGTTGTATAGTCGAGTTTATTTTGTAAAATTAGTTCACGTCTATCAAGAAA'
'''
TCCAATTTTTT
AAATGGAAAACCCATCATCGCACCTGTCGGGCCATCTTCTTTATTAACCAGGTTACTTTGAAGCAATTATTTAATTTTTA
AGAAAAATTAATTCGAAATTTAAGATTACAATCGTTCCTGCCTCAAATGGAAGAAGCCAATCGAAATATGCCGGAAAACG
CGTCAGCTGACGCTTTTCATATTGAAAATGTCGAGGATGACAGCTCGGACGACAGCGATTCGAATTCAGAATCGATTTCG
GAAGAGGTACCTGTCACAACTGACGAAAAGAGCACAGCGAGTAGTACACAGAGGATTGAAATTGATTTAGACGTGTTCAA
GGAGATGAACAGCGCAGTTGATGACAGAGATGTCGGAATTCAGAATGTCGAAGCTCTTCCAGAAGCATTTCAGTCAAAAG
GAGAGGATTCGAAACCGGAAATTACAAAGAAACCGTTGATTGAAGAATTATGAAATTTATGTTTTTTCTATTGTTGTGAA
TTTTTTGTTGAATTTGTTACTTTTTGACCTACCCCTCGAAACATTTTCTCACATTGCGAAATGAATCAATAAAATTTCGA
TTTTAGATTTTTATTCATTTATCTCATTGGTTTCTATCCATTTTTGAATTACCTTCCCATTTTTATTTATCTCTGCCTTT
TCTAAAATTCTCCGAAAACTTTT
'''

from itertools import combinations
import sys
import isoform


#       0        9      16    22       31     38       47  
#                D      A     D        A      A          
seq1 = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
# iso1: AACATGACC        TTACCGTCACATTAGTTCGGAGCCCTATATA
# iso2: AACATGACC                       TTCGGAGCCCTATATA
# iso3: AACATGACC                              CCCTATATA
# iso4: AACATGACCGTTGCGAGTTACC          TTCGGAGCCCTATATA
# iso5: AACATGACCGTTGCGAGTTACC                 CCCTATATA
# iso6: AACATGACC        TTACC          TTCGGAGCCCTATATA
# iso7: AACATGACC        TTACC                 CCCTATATA

# should i count all possible combinations or all possible valid combinations?

seq4 = 'ACACACACGTACACACACACACAGACACACGTACACACCAGACACA'

with open(sys.argv[1], 'r') as ff:

	seq2 = ''
	for line in ff.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			ID = line
		else:
			seq2 += line

# from gff3, to get gene sequence from transcript
# ch.9940	WormBase	gene	102	724	.	+	.	Parent=Transcript:H06H21.11.1
# ch.9940 WormBase        three_prime_UTR
# ch.9940 WormBase        five_prime_UTR

print(seq2[101:110])
print(seq2[495:498])
seq = seq1
print(seq)
# default is 25
minin = 3
minex = 3
minflank5 = 1
minflank3 = 1
maxs = 100
# max number of splices

dons = []
accs = []
for i in range(len(seq)):
	if seq[i:i+2] == 'GT':
		dons.append(i)
	if seq[i:i+2] == 'AG':
		accs.append(i+1)
print(dons, accs)
print('************')
# make an apc algorithm that does not use weights for each individual gene
# make a single set of weights apply to all genes
# how to do?

def short_introns(dons, accs, minin):
	
	for d, a in zip(dons, accs):
		intron_length = a - d + 1
		if intron_length < minin:
			return True
	return False

nshort_introns = 0
nsites = min(len(dons), len(accs), maxs)
for n in range(1, nsites+1):
	for dsites in combinations(dons, n):
		for asites in combinations(accs, n):
			if short_introns(dsites, asites, minin):		
				nshort_introns += 1
				continue
			print(dsites, asites)
			

print('************')

print(len(seq))
dons, accs = (9, 22), (31, 38)

if dons[0] >= minflank5: print('keep')
print(accs[-1])

for d, a in zip(dons, accs):
	if d >= minflank5:
		print('keep')
	

# testing old code
flank = 1
# gets the same sequence each time
#for i in isoform.all_possible(seq,minin,minex,maxs,flank)[0]:
#	print(i)
#print(isoform.all_possible(seq,minin,minex,maxs, flank)[0])

# time with seq3: 0.084s
# WAAAAAAAAAY faster
# i messed up

# using seq1
# introns from old code is the same as new code
# 'introns': [(9, 16)]
# 'introns': [(9, 31)]
# 'introns': [(9, 38)]
# 'introns': [(22, 31)]
# 'introns': [(22, 38)]
# 'introns': [(9, 16), (22, 31)]
# 'introns': [(9, 16), (22, 38)]



























