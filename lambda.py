import math

# need to do a binary search
# look up that algorithm
# for cutoff, determine when number starts changing too much

f = {'A': 0.4, 'C': 0.2, 'G': 0.1, 'T': 0.3}
m = [2, -1]

def qij(l, pi, pj, sij):
	
	return pi * pj * math.exp(l * sij)
	
def qsum(s1, s2, f, m):
	
	for i, j in zip(s1, s2):
		pi = f[i]
		pj = f[j]
		print(pi, pj)
	
print(qij(0.99, 0.2, 0.4, 2))

seq1 = 'AAAACCCCGGGGTTTT'
seq2 = 'ACGTACGTACGTACGT'

qsum(seq1, seq2, f, m)

	
	
