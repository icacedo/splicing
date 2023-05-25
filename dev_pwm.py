import modelib as ml
import sys

seqs = ml.read_txt_seqs(sys.argv[1])

seqs = [
	'ACTCAA',
	'CGTACG',
	'CTAGTG'
]

pfm = [{'A': 0, 'C': 0, 'G': 0, 'T':0} for x in range(len(seqs[0]))]
for i in range(len(seqs)):
	for j in range(len(seqs[i])):
		pfm[j][seqs[i][j]] += 1

print(pfm)
