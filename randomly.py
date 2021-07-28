import isoform
import random

def randseq(n):
	seq = ''
	for i in range(n):
		seq += random.choice('ACGT')
	return seq

min_ex = 15
min_in = 35
trials = 10
minseq = 100
maxseq = 300
step = 50


for n in range(minseq, maxseq+minseq, step):
	for t in range(trials):
		s = randseq(n)
		txs, info = isoform.all_possible(s, min_in, min_ex)
		print(t, n, len(txs), info)