import argparse
import modelib as ml

parser = argparse.ArgumentParser(
	description='generate and score alternative isoforms')
parser.add_argument('fasta', type=str, metavar='<file>',
	help='input single sequence fasta file')
parser.add_argument('--min_intron', required=False, type=int, default=25,
	metavar='<int>', help='minimum length of intron [%(default)i]')
# how does this work? and -help not working


'''
# ch.9940.fa
fp = sys.argv[1]

seqid = None
seq = None
for seqid, seq in ml.read_fastas(fp):
	seqid = seqid
	seq = seq
print(seqid)
print(seq)

dons, accs = ml.get_gtag(seq)

maxs = 100
minin = 25
minex = 25
flank = 100

apc_isoforms, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)
for iso in apc_isoforms:
	print(iso)
'''
