import itertools

seq1 = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
seq3 = 'CTTTAACTGTTTCTCTATTCAGTTGTATAGTCGAGTTTATTTTGTAAAATTAGTTCACGTCTATCAAGAAA'

seq = seq3
minin = 3
minex = 1
maxs = 100
flank = 1

dons = []
accs = []
for i in range(len(seq)):
	if seq[i:i+2] == 'GT':
		dons.append(i)
	if seq[i:i+2] == 'AG':
		accs.append(i+1)
print(dons)

isoforms = []
sites = min(len(dons), len(accs), maxs)
for n in range(1, sites+1):
	for dsites in itertools.combinations(dons, n):
		for asites in itertools.combinations(accs, n):
			print(dsites, asites)
