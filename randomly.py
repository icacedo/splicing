import isoform
from statistics import mean, stdev
import time

trials = 1000
intron = 35
exon = 15
splice = 3
flank = 100

print('\t'.join(('len', 'dons', 'stdev', 'isos', 'stdev', 'time')))
for length in range(300, 1000, 50):
	isos = []
	dons = []
	accs = []
	secs = []
	
	t0 = time.time_ns()
	for t in range(trials):
		s = isoform.randseq(length)
		txs, info = isoform.all_possible(s, intron, exon, splice, flank)
		isos.append(len(txs))
		dons.append(info['donors'])
		accs.append(info['acceptors'])
	t1 = time.time_ns()
	secs = (t1 - t0)/1e9

	print(f'{length}\t{mean(dons):.3f}\t{stdev(dons):.3f}\t{mean(isos):.3f}\t{stdev(isos):.3f}\t{secs:.3f}')
