# triangular smoothing
# new goal: look into curve fitting in python
# fit a gumbel distribution to the length data?
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

m = 5

intsizes = []
for size in intlines:
	intsizes.append(len(size))

intbins = [0 for x in range(max(intsizes))]

for i in range(len(intsizes)):
	intbins[intsizes[i]-1] += 1
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
		#print(bef[j], m2-1+j, '***')
		coefb = m2 - 1 + j
		befx = bef[j] * coefb
		tbefx += befx
		cbefx += coefb
	taftx = 0
	caftx = 0
	for k in range(len(aft)):
		#print(aft[k], m2-k, '@@@')
		coefa = m2 - k
		aftx = aft[k] * coefa
		taftx += aftx
		caftx += coefa
	total_coef = (m2 + 1) + cbefx + caftx
	total = nowx + tbefx + taftx
	smoopt = total/total_coef
	print(bef, now, aft, smoopt)

