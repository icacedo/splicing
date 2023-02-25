#from itertools import combinations

seq1 = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
seq3 = 'CTTTAACTGTTTCTCTATTCAGTTGTATAGTCGAGTTTATTTTGTAAAATTAGTTCACGTCTATCAAGAAA'

seq = seq1
minin = 3

def filter_dup_sites(iso):

	listy = []
	for i in iso:
		for s in i:
			listy.append(s)

	if len(listy) == len(set(listy)):
		return True

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
print('*****')

lengths = []
for i in range(len(introns)):
	lengths.append(len(introns[i]))
	


for i in range(len(introns)):
	for j in range(len(introns[i])):
		iso = []
		for k in range(len(introns)):
			print(i,k,j)
			#introns[k][j]





















'''
isoforms = []
for i in range(len(introns)-1):
	count = 0
	for j in range(len(introns[i])):
		isoforms.append([introns[i][j]])
		isoform = []
		isoform = [None] * len(introns)
		isoform[0] = introns[i][j]
		for intron in introns[i+1]:
			if count == 0:
				isoforms.append([intron])
			isoform[i+1] = intron
			isoforms.append(isoform)

print(isoforms)
# before it was 15
#print(len(isoforms))

#print(set(isoforms[12]))
#print(len(isoforms))
'''



'''
listy = [None] * 4
print(listy)
listii = [[1,2],[3,4]]
print(listii[1])
'''

'''
print('***')
one = ['1', '2', '3']
two = ['11', '22', '33', '44']
three = ['0', '00']
result = zip(one,two,three)
print(list(result))
'''

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

