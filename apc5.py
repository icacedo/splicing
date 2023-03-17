import sys
import itertools

introns1 = [(3,4),(5,6),(7,8)]

'''
def count_down(n):
    print(n)
    if n>0:
        count_down(n-1)

#tator2 = iter([x for x in introns[1:]])
#print([y for y in tator2])

# this works
def afunc(introns):	
	print(introns)
	if len(introns) > 1:
		return afunc(introns[1:])

afunc(introns1)

print('***')
# you don't need a recursive function lol
#for i in range(1, len(introns)):
#	print(introns[i:])

def funk(introns):
	for i in range(1, len(introns)):
		yield introns[i]

for i in funk(introns1):
	print(i)

for i in introns1[0]:
	print(i)
	for j in funk(introns1):
		for k in j:
			print(i,k)

tator = iter([x for x in introns])
print(next(tator))
for i in next(tator):
	print(i)
'''
# next thing to try

# asign each intron in a donor group an acceptor value
# add introns to list only if acceptor value is not already added
# iterate by the longest set of introns by donor
# repeat through shorter intron groups with slicing

introns2 = [[(9, 16), (9, 31), (9, 38)], [(22, 31), (22, 38)]]

introns = introns1

print(introns)
acckey = {}
for i in range(len(introns)):
	for j in introns[i]:
		acckey[j]=i
print(acckey)

listy = []
vals = []
for key in acckey:
	
	if acckey[key] not in vals:
		vals.append(acckey[key])
		listy.append(key)
#print(vals)
print(listy)

acck2 = [(0,1),(2,3)]

def test(acck2, x, y):
	
	for i in range(len(acck2)):	
		print(acck2[i][y])

#test(acck2, 0, 0)
#test(acck2, 0, 1)

ints = [(0,1,2),(3,4),(5,)]

#len(ints[1]) X len(ints[2]) = # of ints[0][0]

times = 1
for i in range(len(ints)-1):
	times *= len(ints[len(ints)-1-i])
	
print(times)

for i in range(len(ints)):
	print(ints[i:])

# print everything but the current index?		
print(ints[1::3])

# to sort into mutuall exclusive lists, have to check if asites
# are before the previous introns acceptor
# can't sort into arrays
#introns2 = [[(9, 16), (9, 31), (9, 38)], [(22, 31), (22, 38)]]
print('******')

dons = [9, 22]
accs = [16, 31, 38]
sites = min(len(dons), len(accs))
intruns = [(9, 16), (9, 31), (9, 38),(22, 31), (22, 38)]
for n in range(1, sites+1):
	for dsites in itertools.combinations(dons, n):
		for asites in itertools.combinations(accs, n):
			print(dsites, asites)

exintrons = [[(9, 16), (22, 31)], [(9, 16), (22, 38)], [(9, 31),], [(9, 38),]]
print('*******')

borb = exintrons

#for x in exintrons:
	#print(x)
	#for i in range(len(x)):
		#for iso in itertools.combinations(x, i):


# still not great	
for x in exintrons:
	for n in range(1, len(x)+1):
		for i in itertools.combinations(x, n):
			print(i)
#[(9, 38), (22, 38)]]
#[(9, 38), (22, 38)]]

print('***********')

introns = [[(9, 16), (9, 31), (9, 38)], [(22, 31), (22, 38)]]
# would still have to check acceptor to donor, too complex
for i in introns:
	dith = {}
	for x in i:
		dith[x]=x[1]

	print(dith)

