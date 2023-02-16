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


#       0        9      16    22       31     38         
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

seq = seq3
print(seq)
# default is 25
minin = 3
# still need to do exons



don = []
acc = []
for i in range(len(seq)):
	if seq[i:i+2] == 'GT':
		don.append(i)
	if seq[i:i+2] == 'AG':
		acc.append(i+1)

introns = []
short_introns = []
for gt in don:
	for ag in acc:
		if gt > ag:
			continue
		if ag - gt + 1 < minin:
			short_introns =+ 1
			continue
		introns.append((gt,ag))
print(introns)

print('************')

# to speed up, maybe don't generate sites in the APC loop?
# try using previously coded sanity checks instead of building them in?
isoforms = []
for i in range(1, len(introns)+1):
	for com in combinations(introns, i):
		sites = []
		for k,j in com:
			sites += k,j
		if len(sites) != len(set(sites)):
			continue
		elif sites != sorted(sites):
			continue
		else:
			isoforms.append(com)
			
# need to get isoform information, where exons and introns and don/acc sites are
			

print(isoforms)

# time with seq3: real 28.582s
# this code takes forever
# need to see how long it took the old code run

'''
print('************')
# testing old code
minex = 1
maxs = 100
flank = 1
# gets the same sequence each time
#for i isoform.all_possible(seq,minin,minex,maxs,flank)[0]:
#	print(i)
print(isoform.all_possible(seq,minin,minex,maxs, flank)[0])
'''
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



























