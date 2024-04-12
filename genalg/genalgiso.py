import argparse
import random
from datetime import datetime
import apc_model_lib as aml
import os

parser = argparse.ArgumentParser(
	description='genetic algorithm for weight and icost optimization')
parser.add_argument('fasta', type=str, metavar='<file>',
	help='single gene fasta file')
parser.add_argument('gff', type=str, metavar='<file>',
	help='single gene gff file')
parser.add_argument('--program', required=False, type=str, 
	default='../apc/apc_isogen.py', metavar='<exec>', 
	help='path to apc program [%(default)s]')

parser.add_argument('--max_splice', required=False, type=int, default=3,
	metavar='<int>', help='maximum number of splicing events %(default)d')
parser.add_argument('--min_intron', required=False, type=int, default=25,
	metavar='<int>', help='minimum length of intron %(default)d')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length of exon %(default)d')
parser.add_argument('--flank', required=False, type=int, default=100,
	metavar='<int>', help='length of genomic flank on each side %(default)d')
parser.add_argument('--limit', required=False, type=int, default=100, 
	metavar='<int>', help='limit number of saved apc isoforms %(default)d')

parser.add_argument('--elen', required=False, type=str, 
	default='../mkmdls_out/exon_len.tsv', metavar='<file>', 
	help='path to exon length model .tsv [%(default)s]')
parser.add_argument('--ilen', required=False, type=str, 
	default='../mkmdls_out/intron_len.tsv', metavar='<file>',
	help='path to intron length model .tsv [%(default)s]')
parser.add_argument('--emm', required=False, type=str, 
	default='../mkmdls_out/exon_mm.tsv', metavar='<file>',
	help='path to exon markov model .tsv [%(default)s]')
parser.add_argument('--imm', required=False, type=str, 
	default='../mkmdls_out/intron_mm.tsv', metavar='<file>',
	help='path to intron markov model .tsv [%(default)s]')
parser.add_argument('--dpwm', required=False, type=str, 
	default='../mkmdls_out/donor_pwm.tsv', metavar='<file>',
	help='path to donor pwm .tsv [%(default)s]')
parser.add_argument('--apwm', required=False, type=str, 
	default='../mkmdls_out/acceptor_pwm.tsv', metavar='<file>',
	help='path to acceptor pwm .tsv [%(default)s]')
parser.add_argument('--icost', required=False, type=float, default=0.0,
	metavar='<float>', help='intron cost %(default).2d')
	
parser.add_argument('--pop', required=False, type=int, default=50,
	metavar='<int>', help='population size [%(default)i]')
parser.add_argument('--gen', required=False, type=int, default=50,
	metavar='<int>', help='number of generations [%(default)i]')
parser.add_argument('--die', required=False, type=float, default=0.5,
	metavar='<float>', help='fraction tha die each generation [%(default).2f]')

args = parser.parse_args()

random.seed(datetime.now().timestamp())

def chrom():

	genotype = {
		'--welen': random.random(),
		'--wilen': random.random(),
		'--wemm': random.random(),
		'--wimm': random.random(),
		'--wdpwm': random.random(),
		'--wapwm': random.random(),
		'--icost': random.randint(0, 100)
	}

	return genotype

def bad_introns(gff):

	with open(gff, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			if line[1] == 'RNASeq_splice': 
				beg = int(line[3])
				if beg <= 100 + args.min_exon: continue
				print(line)

bad_introns(args.gff)




'''
def get_fit(chrom):

	tmpfile = f'tmp.{os.getpid()}.gff'
	cmd = (
		f'python3 {args.program} {args.fasta} --gff {args.gff} '
		f'--max_splice {args.max_splice} --min_intron {args.min_intron} '
		f'--min_exon {args.min_exon} --flank {args.flank} ' 
		f'--limit {args.limit} --exon_len {args.elen} '
		f'--intron_len {args.ilen} --exon_mm {args.emm} '
		f'--intron_mm {args.imm} --donor_pwm {args.dpwm} '
		f'--acceptor_pwm {args.apwm} > {tmpfile}'
		)

	os.system(cmd)

	with open(tmpfile, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			

get_fit(chrom)
'''


















