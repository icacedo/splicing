import argparse
import modelib as ml
import csv
import sys
import math
import json

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

exon_scores = {}
intron_scores = {}
for iso in apc_isoforms:
	total_iso_score = 0
	for exon in iso['exons']:	
		if exon in exon_scores: continue
		elen_score = ml.get_exin_len_score(exon, re_elen_log2, ea, eb, eg)
		emm_score = ml.get_exin_mm_score(exon, seq, re_emm_log2)
		escore = elen_score + emm_score
		exon_scores[exon] = escore
	for intron in iso['introns']:
		if intron in intron_scores: continue
		ilen_score = ml.get_exin_len_score(intron, re_ilen_log2, ia, ib, ig)
		imm_score = ml.get_exin_mm_score(intron, seq, re_imm_log2, 'GT', 'AG')
		dseq, aseq = ml.get_donacc_seq(intron, seq)
		dpwm_score = ml.get_donacc_pwm_score(dseq, re_dpwm)
		apwm_score = ml.get_donacc_pwm_score(aseq, re_apwm)
		iscore = ilen_score + imm_score + dpwm_score + apwm_score
		intron_scores[intron] = iscore
	for exon in iso['exons']:
		total_iso_score += exon_scores[exon]
	for intron in iso['introns']:
		total_iso_score += intron_scores[intron]
	total_iso_score -= len(iso['introns']) * args.icost
	iso['score'] = total_iso_score

apc_isoforms = sorted(apc_isoforms, key=lambda iso: iso['score'], reverse=True)

iso_weights = []
iso_total = 0
for iso in apc_isoforms:
	iso_weight = 2 ** iso['score']
	iso_weights.append(iso_weight)
	iso_total += iso_weight

iso_probs = []
for w in iso_weights:
	iso_probs.append(w / iso_total)

exon_counts = {}
intron_counts = {}
exon_total = 0
intron_total = 0
for iso in apc_isoforms:
	for exon in iso['exons']:
		if exon not in exon_counts:
			exon_counts[exon] = 1
			exon_total += 1
		else:
			exon_counts[exon] += 1
			exon_total += 1
	for intron in iso['introns']:
		if intron not in intron_counts:
			intron_counts[intron] = 1
			intron_total += 1
		else:
			intron_counts[intron] += 1
			intron_total += 1

exon_freqs = {}
intron_freqs = {}
for exon in exon_counts:
	exon_freqs[exon] = exon_counts[exon] / exon_total
for intron in intron_counts:
	intron_freqs[intron] = intron_counts[intron] / intron_total

name = seqid.split(' ')[0]

print('# name:', name)
print('# length:', len(seq))
print('# donors:', len(dons))
print('# acceptors:', len(accs))
print('# icost:', args.icost)
print('# trials:', trials)
print('# isoforms:', len(apc_isoforms))
print('# complexity:', f'{ml.get_entropy(iso_probs):.4f}')

gff_writer = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
gff_writer.writerow([name, 'apc_isogen', 'gene', iso['beg']+1, iso['end']+1,
	'.', '+', '.', 'ID=Gene-' + name])
gff_writer.writerow([])
count = 0
for iso in apc_isoforms:
	if count <= args.limit - 1:
		iso_prob_f = '{:.5e}'.format(iso_probs[count])
		gff_writer.writerow([name, 'apc_isogen', 'mRNA', iso['beg']+1, 
			iso['end']+1, iso_prob_f, '+', '.', 'ID=iso-'+name+'-'+
			str(count+1)+';Parent=Gene-'+name])
		for exon in iso['exons']:
			escore_f = '{:.5e}'.format(exon_scores[exon])
			efreq_f = '{:.5e}'.format(exon_freqs[exon])
			gff_writer.writerow([name, 'apc_isogen', 'exon', exon[0]+1,
				exon[1]+1, iso_prob_f, '+', '.', 'Parent='+'iso-'+name+'-'
				+str(count+1)+';score='+str(escore_f)+';exfreq='+str(efreq_f)])
		for intron in iso['introns']:
			iscore_f = '{:.5e}'.format(intron_scores[intron])
			ifreq_f = '{:.5e}'.format(intron_freqs[intron])
			gff_writer.writerow([name, 'apc_isogen', 'intron', intron[0]+1,
				intron[1]+1, iso_prob_f, '+', '.', 'Parent='+'iso-'+name+'-'
				+str(count+1)+';score='+str(iscore_f)+';infreq='+str(ifreq_f)])
		gff_writer.writerow([])
		count += 1



