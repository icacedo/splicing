from itertools import combinations
import sys
import isoform

'''
with open(sys.argv[1], 'r') as ff:

	seq = ''
	for line in ff.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			ID = line
		else:
			seq += line

# minimum intorn/exon size, default is 25
minin = 3
minex = 4
# genomic flank, default 100
flank = 5
# max number of splices, default 100
maxs = 100
'''
##### Some notes #####
'''
make an apc algorithm that does not use weights for each individual gene
make a single set of weights apply to all genes
how to do?
not sure i see the utility in counting discarded isoforms
by short exons on introns
if discarded due to short exon first
but there is a short intron
then it will only count it as a short exon?
maybe just go by discarded isoforms
'''

def get_gtag(seq):

	dons = []
	accs = []
	for i in range(len(seq)):
		if seq[i:i+2] == 'GT':
			dons.append(i)
		if seq[i:i+2] == 'AG':
			accs.append(i+1)

	return dons, accs

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

def apc(dons, accs, maxs, minin, minex, flank):
	
	apc_isoforms = []

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
				if short_introns(dsites, asites, minin): continue
				if short_exons(dsites, asites, flank, minex): continue
				apc_isoform['seq'] = seq
				apc_isoform['beg'] = flank
				apc_isoform['end'] = len(seq) - flank - 1
				apc_isoform['exons'] = get_exons(dsites, asites, flank, seq)
				apc_isoform['introns'] = get_introns(dsites, asites)
				print(apc_isoform)
	print(trials)














