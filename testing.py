from itertools import combinations

#       0        9      16    22       31     38       47  
#                D      A     D        A      A          
seq = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
# iso1: AACATGACC        TTACCGTCACATTAGTTCGGAGCCCTATATA
# iso2: AACATGACC                       TTCGGAGCCCTATATA
# iso3: AACATGACC                              CCCTATATA
# iso4: AACATGACCGTTGCGAGTTACC          TTCGGAGCCCTATATA
# iso5: AACATGACCGTTGCGAGTTACC                 CCCTATATA
# iso6: AACATGACC        TTACC          TTCGGAGCCCTATATA
# iso7: AACATGACC        TTACC                 CCCTATATAj

minin = 3
minex = 5
flank = 4
maxs = 100

def get_gtag(seq):

	dons = []
	accs =[]
	for i in range(len(seq)):
		if seq[i:i+2] == 'GT':
			dons.append(i)
		if seq[i:i+2] == 'AG':
			accs.append(i+1)
		
	return dons, accs

def short_introns(dons, accs, minin):
	
	for d, a in zip(dons, accs):
		intron_length = a - d + 1
		if intron_length < minin:
			return True
		
	return False

# (9, 22) (16, 38)
def short_exons(dons, accs, flank, minex):
		
	fexlen = dons[0] - flank
	if fexlen < minex:
		return True

	lexbeg = accs[-1] + 1
	lexend = len(seq) - flank
	lexlen = lexend - lexbeg + 1
	if lexlen < minex:
		return True

	for i in range(len(dons)-1):
		if dons[i+1] - asites[i] < minex:
			return True

	return False



dons, accs = get_gtag(seq)
print(dons, accs)
print('**********')

nsites = min(len(dons), len(accs), maxs)
for n in range(1, nsites+1):
	for dsites in combinations(dons, n):
		for asites in combinations(accs, n):
			if short_introns(dsites, asites, minin):
				continue
			if short_exons(dsites, asites, flank, minex):
				continue
			print(dsites, asites)
			#print(dsites, asites, '***')	
			#print(short_introns(dsites, asites, minin), 'intron')
			#print(short_exons(dsites, asites, flank, minex), 'exon')


