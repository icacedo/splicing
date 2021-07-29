import isoform
import random
import argparse
from statistics import mean, stdev
import time

def randseq(n):
	seq = ''
	for i in range(n):
		seq += random.choice('ACGT')
	return seq

if __name__ == '__main__':

	parser = argparse.ArgumentParser(
		description='Random experiments with isoforms')
	parser.add_argument('length', type=int, metavar='<length>',
		help='length of random sequence')                
	parser.add_argument('--intron', required=False, type=int, default=35,
		metavar='<int>', help='minimum length of intron [%(default)i]')
	parser.add_argument('--exon', required=False, type=int, default=25,
		metavar='<int>', help='minimum length exon [%(default)i]')
	parser.add_argument('--trials', required=False, type=int, default=100,
		metavar='<int>', help='number of runs [%(default)i]')
	arg = parser.parse_args()

isos = []
dons = []
accs = []
secs = []
for t in range(arg.trials):
	s = randseq(arg.length)
	t0 = time.time_ns()
	txs, info = isoform.all_possible(s, arg.intron, arg.exon)
	t1 = time.time_ns()
	isos.append(len(txs))
	dons.append(info['donors'])
	accs.append(info['acceptors'])
	secs.append((t1 - t0)/1e9)

#print(f'isos {mean(isos):.3f} {stdev(isos):.3f} {min(isos)} {max(isos)}')
#print(f'dons {mean(dons):.3f} {stdev(dons):.3f} {min(dons)} {max(dons)}')
#print(f'accs {mean(accs):.3f} {stdev(accs):.3f} {min(accs)} {max(accs)}')
#print(f'time {mean(secs):.3f} {stdev(secs):.3f} {min(secs)} {max(secs)}')

print(f'{arg.length}\t{stdev(dons):.3f}\t{stdev(accs):.3f}\t{mean(isos):.3f}\t{mean(secs):.3f}')
