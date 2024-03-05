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
'''
with open(sys.argv[1], 'r') as ff:

	seq2 = ''
	for line in ff.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			ID = line
		else:
			seq2 += line
'''
# from gff3, to get gene sequence from transcript
# ch.9940	WormBase	gene	102	724	.	+	.	Parent=Transcript:H06H21.11.1
# ch.9940 WormBase        three_prime_UTR
# ch.9940 WormBase        five_prime_UTR

#print(seq2[101:110])
#print(seq2[495:498])
seq5 = 'CCCCCAAAAAGTCCCCCCAGAAAAACCCCC'
# 5 nt flank, 5nt exon, 10nt intron, 5nt exon, 5nt flank
# 0-4 flank, 5-9 exon, 10-19 intron, 20-24 exon, 25-29 flank
seq = seq5
print(seq)
# default is 25
minin = 3
minex = 4
# default should be 100
flank = 5
maxs = 100
# max number of splices

def get_gtag(seq):

	dons = []
	accs = []
	for i in range(len(seq)):
		if seq[i:i+2] == 'GT':
			dons.append(i)
		if seq[i:i+2] == 'AG':
			accs.append(i+1)

	return dons, accs

dons, accs = get_gtag(seq)
print('************')
# make an apc algorithm that does not use weights for each individual gene
# make a single set of weights apply to all genes
# how to do?

# using index starting at 0
def short_introns(dons, accs, minin):
	
	for d, a in zip(dons, accs):
		intron_length = a - d + 1
		if intron_length < minin:
			return True

	return False

# using index starting at 0
def short_exons(dons, accs, flank, minex):
	
	# check 5' first exon
	fexlen = dons[0] - flank
	if fexlen < minex:
		return True
	
	# check 3' last exon
	lexbeg = accs[-1] + 1
	lexend = len(seq) - flank - 1 
	lexlen = lexend - lexbeg + 1
	if lexlen < minex:
		return True
	
	# check interior exons
	for i in range(len(dons)-1):
		if dons[i+1] - accs[i] - 1 < minex:
			return True

	return False	
	
def get_exons(dsites, asites, flank, seq):

	exons = []
	exons.append((flank, dsites[0]-1))
	for i in range(1, len(dsites)):
		exbeg = asites[i-1] + 1
		exend = dsites[i] - 1
		exons.append((exbeg, exend))
	exons.append((asites[-1]+1, len(seq)-flank-1))

	return exons	

def get_introns(dsites, asites):

	introns = []	
	for d, a in zip(dsites, asites):	
		introns.append((d, a))
	
	return introns

# not sure i see the utility in counting discarded isoforms
# by short exons on introns
# if discarded due to short exon first
# but there is a short intron
# then it will only count it as a short exon?
# maybe just go by discarded isoforms


apc_isoform = {
	'seq': '',
	'beg': '',
	'end': '',
	'exons': [],
	'introns': [],
	'score': 0
}

trials = 0
short_introns_exons = 0
nsites = min(len(dons), len(accs), maxs)
for n in range(1, nsites+1):
	for dsites in combinations(dons, n):
		for asites in combinations(accs, n):
			trials += 1
			if short_introns(dsites, asites, minin):
				continue
			if short_exons(dsites, asites, flank, minex):
				continue
			apc_isoform['seq'] = seq
			apc_isoform['beg'] = flank
			apc_isoform['end'] = len(seq) - flank - 1
			apc_isoform['exons'] = get_exons(dsites, asites, flank, seq)
			apc_isoform['introns'] = get_introns(dsites, asites)
			print(apc_isoform)
print(trials)
print('************')
# testing old code
# old code offset by 1 at the end, i believe this is a mistake
#flank = 1
# gets the same sequence each time
for i in isoform.all_possible(seq,minin,minex,maxs,flank)[0]:
	print(i)
print(isoform.all_possible(seq,minin,minex,maxs, flank)[0])

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




















