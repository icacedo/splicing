import argparse
import modelib as ml
import csv
import sys
import math

# params for test seq: maxs 100, minin 3, minex 4, flank 5
parser = argparse.ArgumentParser(
	description='generate and score alternative isoforms')
parser.add_argument('fasta', type=str, metavar='<file>',
	help='input single sequence fasta file')
parser.add_argument('--gff', type=str, metavar='<file>', required=False,
	help='input .gff3 for single gene')

# apc parameters
parser.add_argument('--max_splice', required=False, type=int, default=3,
	metavar='<int>', help='maximum number of splicing events %(default)d')
parser.add_argument('--min_intron', required=False, type=int, default=25,
	metavar='<int>', help='minimum length of intron %(default)d')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length of exon %(default)d')
parser.add_argument('--flank', required=False, type=int, default=100,
	metavar='<int>', help='length of genomic flank on each side %(default)d')
parser.add_argument('--limit', required=False, type=int, default=20, 
	metavar='<int>', help='limit number of saved apc isoforms %(default)d')

# probabilistic models
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

seqid = None
seq = None
for seqid, seq in ml.read_fastas(args.fasta):
	seqid = seqid
	seq = seq

if args.gff:
	dons, accs = ml.read_gff_sites(seq, args.gff) 
else:
	dons, accs = ml.get_gtag(seq)

maxs = args.max_splice
minin = args.min_intron
minex = args.min_exon
flank = args.flank

apc_isoforms, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)

if args.exon_len:
	re_elen_pdf, re_elen_log2 = ml.read_exin_len(args.exon_len)
	ea, eb, eg = ml.read_len_params(args.exon_len) 
if args.intron_len:
	re_ilen_pdf, re_ilen_log2 = ml.read_exin_len(args.intron_len)
	ia, ib, ig = ml.read_len_params(args.intron_len)
if args.exon_mm:
	re_emm_prob, re_emm_log2 = ml.read_exin_mm(args.exon_mm)
if args.intron_mm:
	re_imm_prob, re_imm_log2 = ml.read_exin_mm(args.intron_mm)
if args.donor_pwm:
	re_dppm, re_dpwm = ml.read_pwm(args.donor_pwm)
if args.acceptor_pwm:
	re_appm, re_apwm = ml.read_pwm(args.acceptor_pwm)

def get_exin_len(exin):

	exin_len = exin[1] - exin[0] + 1
	return exin_len	

def get_exin_len_score(exin, exin_len_model, a, b, g):

	exin_len = get_exin_len(exin)

	if exin_len < len(exin_len_model):
		exin_score = exin_len_model[exin_len]
	else:
		exin_prob = ml.frechet_pdf(length, a, b, g)
		expect = 1/exin_len
		exin_score = math.log2(exin_prob/expect)
	return exin_score

def get_exin_seq(exin, seq):

	beg = exin[0]
	end = exin[1] + 1
	exin_seq = seq[beg:end]
	return exin_seq

def get_exin_mm_score(exin, seq, exin_mm, dpwm=None, apwm=None):

	exin_seq = get_exin_seq(exin, seq)

	k = 0
	for key in exin_mm:
		k = len(key)
		break

	if dpwm and apwm:
		exin_seq = exin_seq[len(dpwm):-len(apwm)]
	
	exin_score = 0
	for i in range(len(exin_seq)):
		if len(exin_seq[i:i+k]) == k:
			kmer = exin_seq[i:i+k]
			exin_score += float(exin_mm[kmer])

	return exin_score
	

exons = {}
for iso in apc_isoforms:
	for exon in iso['exons']:
		len_score = get_exin_len_score(exon, re_elen_log2, ea, eb, eg)
		mm_score = get_exin_mm_score(exon, seq, re_emm_log2)
		print(len_score, mm_score)
		
	




