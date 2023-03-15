introns = [(3,4),(5,6),(7,8)]

def count_down(n):
    print(n)
    if n>0:
        count_down(n-1)

count_down(5)


tator2 = iter([x for x in introns[1:]])
print([y for y in tator2])

# this works
def afunc(introns):
	print(introns)
	if len(introns)>1:
		afunc(introns[1:])

afunc(introns)


'''
def afunc(introns):
	print(introns)	
	if len(introns) > 1:
		tator = iter([x for x in introns[1:]])
		afunc(tator)
			
print(afunc(introns))
'''
