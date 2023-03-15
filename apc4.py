from itertools import combinations

seq1 = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
seq3 = 'CTTTAACTGTTTCTCTATTCAGTTGTATAGTCGAGTTTATTTTGTAAAATTAGTTCACGTCTATCAAGAAA'

seq = seq1
minin = 3
minex = 1
maxs = 100
flank = 1

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

set1 = (0,1,2)
set2 = (3,4)
set3 = (5,6,7)
set4 = [set1, set2]
set5 = [set1, set2, set3]
'''
def nexter(group):
	test = [x for x in group]
	tator = iter(test)
	for i in range(len(test)):
		yield next(tator)
'''
pik = set4

for i in range(len(pik)):
	#tator = iter([x for x in pik])
	count = 0
	for j in pik[i]:
		tator = iter([x for x in pik])
		#print(j, next(tator))
		if count == 0:
			next(tator)
			count += 1
		for k in next(tator):	
			print(j,k)
print('*****')

for i in range(len(pik)):
	for j in pik[i]:
		tator = iter([x for x in pik[1:]])
		#print(j, next(tator))
		for k in next(tator):
			print(j,k)
print('******')
# this what i was trying to get, but still need to add to a list
for i in pik[0]:
	tator2 = iter([x for x in pik[1:]])
	for j in next(tator2):
		print(i,j)



	
