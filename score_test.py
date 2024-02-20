import modelib as ml
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--exon_len', required=False, type=str, metavar='<file>', 
	help='exon length model .tsv')
parser.add_argument('--intron_len', required=False, type=str, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--exon_mm', required=False, type=str, metavar='<file>',
	help='exon markov model .tsv')
parser.add_argument('--intron_mm', required=False, type=str, metavar='<file>',
	help='intron markov model .tsv')
parser.add_argument('--donor_pwm', required=False, type=str, metavar='<file>',
	help='donor pwm .tsv')
parser.add_argument('--acceptor_pwm', required=False, type=str, metavar='<file>',
	help='acceptor pwm .tsv')
parser.add_argument('--icost', required=False, type=float, default=0.0,
	metavar='<float>', help='intron cost %(default).2d')
	
args = parser.parse_args()

'''
fasta = sys.argv[1]

seq = ''
with open(fasta, 'r') as fp:
	for line in fp.readlines():
		seq += line.rstrip()
'''

seq = 'AGCGAAGTACAAGCAGCGGTAAAGACTAAATGAGTTCTACAAGTTCGAGTAAGATAATACACA'

dons, accs = ml.get_gtag(seq)

maxs = 3
minin = 5
minex = 5
flank = 10

apc_isoforms, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)

print(apc_isoforms)
