import math

# need to do a binary search
# look up that algorithm
# for cutoff, determine when number starts changing too much

freq = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
match = 1
mismatch = -1
 
def qij(l, pi, pj, sij):
	
	return pi * pj * math.exp(l * sij)
	
high = 2
low = 0

while high - low > 1e-6:
	l = (high + low)/2
	print(l)
	qsum = 0
	for nt1, p1 in freq.items():
		for nt2, p2 in freq.items():
			if nt1 == nt2:
				qsum += qij(l, p1, p2, match)
			else:
				qsum += qij(l, p1, p2, mismatch)
	if qsum > 1: 
		high = l
	else:
		low = l

		 

	
	
