# triangular smoothing
# new goal: look into curve fitting in python
# fit a gumbel distribution to the length data?

import argparse
import modelib as ml
import sys


'''
parser = argparse.ArgumentParser()
parser.add_argument('--out_file', type=str, required=False)
arg = parser.parse_args()
'''

intlines = [
	'AAATGCTATA',
	'ATCTA',
	'CGACTCT',
	'ACTACGTACG',
	'ACTAGG',
	'CGCGATCGTT',
	'AGCGATTTACAGCG',
	'AGCAC',
	'CTACTG',
	'AGAG',
	'AGCT'
]

fp = sys.argv[1]

m = 5

intbins = ml.get_intbins(fp)[0]

print(intbins)

m2 = int((m/2) + 0.5 - 1)
for i in range(len(intbins)):
	inx = i-m2
	if inx < 0:
		inx = 0
	bef = intbins[inx:i]
	now = intbins[i]
	aft = intbins[i+1:i+m2+1]
	nowx = now * (m2 + 1)
	tbefx = 0
	cbefx = 0
	for j in range(len(bef)):
		coefb = m2 - 1 + j
		befx = bef[j] * coefb
		tbefx += befx
		cbefx += coefb
	taftx = 0
	caftx = 0
	for k in range(len(aft)):
		coefa = m2 - k
		aftx = aft[k] * coefa
		taftx += aftx
		caftx += coefa
	total_coef = (m2 + 1) + cbefx + caftx
	total = nowx + tbefx + taftx
	smoopt = total/total_coef
	print(bef, now, aft, total, smoopt)
	# index is wrong
	#print('{0},{1}'.format(i, smoopt))


'''
if arg.out_file is not None:
	with open('arg.out_file', 'w') as csvfile:
		csvwriter = csv.writer(csvfile:
		for row in result.items():
'''
