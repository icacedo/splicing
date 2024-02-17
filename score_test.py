import modelib as ml
import sys

fasta = sys.argv[1]

seq = ''
with open(fasta, 'r') as fp:
	for line in fp.readlines():
		seq += line.rstrip()

dons, accs = ml.get_gtag(seq)

maxs = 3
minin = 25
minex = 25
flank = 100

apc_isoforms, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)

print(apc_isoforms)
