#from itertools import combinations

seq1 = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
seq3 = 'CTTTAACTGTTTCTCTATTCAGTTGTATAGTCGAGTTTATTTTGTAAAATTAGTTCACGTCTATCAAGAAA'

seq = seq1
minin = 3

dons = []
accs = []
for i in range(len(seq)):
	if seq[i:i+2] == 'GT':
		dons.append(i)
	if seq[i:i+2] == 'AG':
		accs.append(i+1)

introns = []
short_introns = []
for gt in dons:
	intdons = []
	for ag in accs:
		if gt > ag:
			continue
		if ag - gt + 1 < minin:
			short_introns =+ 1
			continue
		intdons.append((gt,ag))
	introns.append(intdons)
print(introns)

for i in introns:
	dons = iter(i)
	for j in range(len(i)):
		print(next(dons))
	print('***')



'''
isoforms = []
for i in range(1, len(dons)+1):
	for com in combinations(introns, i):
		print(com)
		sites = []
		for k,j in com:
			sites += k,j
		if len(sites) != len(set(sites)):
			continue
		elif sites != sorted(sites):
			continue
		else:
			isoforms.append(com)
'''

#print(isoforms)







# correct output for seq1:
# [((9, 16),), ((9, 31),), ((9, 38),), ((22, 31),), ((22, 38),), ((9, 16), (22, 31)), ((9, 16), (22, 38))]

