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
set3 = (5,6)
set4 = (7,8)
set5 = (9,10,11)
set6 = [set1, set2]
set7 = [set1, set2, set3]
set8 = [set2, set3, set4]
'''
def nexter(group):
	test = [x for x in group]
	tator = iter(test)
	for i in range(len(test)):
		yield next(tator)
'''
pik = set8

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
print('*********')
# this is the best one yet
tator3 = iter([x for x in pik[1:]])

for i in range(len(pik)-1):
	store = next(tator3)
	dct = {}
	for j in pik[0]:
		k_list = []
		for k in store:
			print(j,k)
			k_list.append(k)
			dct[j]=k_list
	print(dct)
print('******')
# when using pik = set8
dct_list = [{3: [5, 6], 4: [5, 6]}, {3: [7, 8], 4: [7, 8]}]

for i in range(len(dct_list)):
	dctlist3 = []
	for j in dct_list[i]:
		for k in dct_list[i][j]:
			newdct = {}
			newdct[j]=k
			print(newdct)
			

'''
for i in dct_list:
	print(i)
	for j in i:
		listy = []
		listy.append(j)
		print(i[j])
		for k in i[j]:
			print(k)
			listy.append(k)
			print(listy, '*')
print('*********')

#tator4 = iter([x for x in dct_list])
#print(next(tator4))
dct2 = {}
for i in range(len(dct_list)):
	for j in dct_list[i]:
		print(j, '*')
		tator4 = iter([x for x in dct_list])
		for k in range(len(dct_list)):
			become = next(tator4)
			for l in become[j]:
				print(j,l)
				#dct2[l]
#print(dct2)	
'''	
print('******')
# each number represents and individual intron
# group by the same donor site
fkint = [[(3), (4)], [(5), (6)], [(7), (8)]]
# maybe represent each intron as a number? then sort

	




	



