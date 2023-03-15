import sys

introns = [(3,4),(5,6),(7,8)]

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

afunc(introns)

print('***')
# you don't need a recursive function lol
#for i in range(1, len(introns)):
#	print(introns[i:])

def funk(introns):
	for i in range(1, len(introns)):
		yield introns[i]

for i in funk(introns):
	print(i)

for i in introns[0]:
	print(i)
	for j in funk(introns):
		for k in j:
			print(i,k)

tator = iter([x for x in introns])
print(next(tator))
for i in next(tator):
	print(i)



