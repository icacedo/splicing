import modelib as ml
import sys
import math 

seqs = ml.read_txt_seqs(sys.argv[1])
'''
seqs = [
	'ACTCAA',
	'CGTACG',
	'CTAGTG',
	'GTCGTA'
]
'''
pfm = [{'A': 0, 'C': 0, 'G': 0, 'T':0} for x in range(len(seqs[0]))]
for i in range(len(seqs)):
	for j in range(len(seqs[i])):
		pfm[j][seqs[i][j]] += 1
print(pfm)
print('*****')
ppm = [{'A': 0, 'C': 0, 'G': 0, 'T':0} for x in range(len(pfm))]
for i in range(len(pfm)):
	for n in pfm[i]:
		ppm[i][n] = pfm[i][n]/len(seqs)
print(ppm)
print('*****')
pwm = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(len(ppm))]
for i in range(len(ppm)):	
	for n in ppm[i]:
		if ppm [i][n] == 0:
			pwm[i][n] = -100
		else:
			pwm[i][n] = math.log2(ppm[i][n]/0.25)
print(pwm)
		
