import argparse
import isomod as im
import csv
import sys
import math
import json
import itertools

parser = argparse.ArgumentParser(
	description='generate and score alternative isoforms')
parser.add_argument('fasta', type=str, metavar='<file>',
	help='input single sequence fasta file')
parser.add_argument('--gff', type=str, metavar='<file>', required=False,
	help='input .gff3 for single gene')

# apc parameters
parser.add_argument('--maxs', required=False, type=int, default=3,
	metavar='<int>', help='maximum number of splicing events %(default)d')
parser.add_argument('--minin', required=False, type=int, default=35,
	metavar='<int>', help='minimum length of intron %(default)d')
parser.add_argument('--minex', required=False, type=int, default=25,
	metavar='<int>', help='minimum length of exon %(default)d')
parser.add_argument('--flank', required=False, type=int, default=99,
	metavar='<int>', help='length of genomic flank on each side %(default)d')
parser.add_argument('--limit', required=False, type=int, default=20, 
	metavar='<int>', help='limit number of saved apc isoforms %(default)d')

# probabilistic models
parser.add_argument('--elen', required=False, type=str, metavar='<file>', 
	help='exon length model .tsv')
parser.add_argument('--ilen', required=False, type=str, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--emm', required=False, type=str, metavar='<file>',
	help='exon markov model .tsv')
parser.add_argument('--imm', required=False, type=str, metavar='<file>',
	help='intron markov model .tsv')
parser.add_argument('--dpwm', required=False, type=str, metavar='<file>',
	help='donor pwm .tsv')
parser.add_argument('--apwm', required=False, type=str, metavar='<file>',
	help='acceptor pwm .tsv')

# penalties
parser.add_argument('--welen', required=False, type=float, metavar='<float>', 
	default=1.0, help='exon length model weight [%(default).2f]')
parser.add_argument('--wilen', required=False, type=float, metavar='<float>', 
	default=1.0, help='intron length model weight [%(default).2f]')
parser.add_argument('--wemm', required=False, type=float, metavar='<float>', 
	default=1.0, help='exon Markov model weight [%(default).2f]')
parser.add_argument('--wimm', required=False, type=float, metavar='<float>', 
	default=1.0, help='intron Markov model weight [%(default).2f]')
parser.add_argument('--wdpwm', required=False, type=float, metavar='<float>', 
	default=1.0, help='donor pwm model weight [%(default).2f]')
parser.add_argument('--wapwm', required=False, type=float, metavar='<float>', 
	default=1.0, help='acceptor pwm model weight [%(default).2f]')
parser.add_argument('--icost', required=False, type=float, default=0.0,
	metavar='<float>', help='intron cost %(default).2d')
	
args = parser.parse_args()
# 13301 is the shortest gene, quick to test
seqid, seq = im.read_fasta(args.fasta)
seq_info = seqid.split(' ')
coor = seq_info[1]
strand = seq_info[2]
wbgene = seq_info[3].split(':')[1]

if args.gff:
	dons, accs = im.read_gff_sites(seq, args.gff) 
else:
	dons, accs = im.get_gtag(seq, args.flank, args.minex)

abc_isoforms, trials = im.abc(dons, accs, args.maxs, args.minin, 
							  args.minex, args.flank, seq)

re_elen = im.read_len(args.elen) if args.elen else None
re_ilen = im.read_len(args.ilen) if args.ilen else None
re_emm = im.read_mm(args.emm) if args.emm else None
re_imm = im.read_mm(args.imm) if args.imm else None
re_dpwm = im.read_pwm(args.dpwm) if args.dpwm else None
re_apwm = im.read_pwm(args.apwm) if args.apwm else None

escores = {}
iscores = {}
dscores = {}
ascores = {}
for iso in abc_isoforms:
	for exon in iso['exons']:
		if exon in escores: continue
		if args.elen: 
			elen_score = im.score_len(re_elen, exon) * args.welen
		else:
			elen_score = 0
		if args.emm: 
			emm_score = im.score_mm(re_emm, exon, seq) * args.wemm
		else:
			emm_score = 0
		escores[exon] = elen_score + emm_score
	for intron in iso['introns']:
		if intron in iscores: continue
		if args.ilen: 
			ilen_score = im.score_len(re_ilen, intron) * args.wilen
		else:
			ilen_score = 0
		if args.imm:
			imm_score = im.score_mm(re_imm, intron, seq, re_dpwm, re_apwm) \
				* args.wimm
		else:
			imm_score = 0
		dseq, aseq = im.get_daseq(intron, seq)
		if args.dpwm: 
			dpwm_score = im.score_pwm(dseq, re_dpwm) * args.wdpwm
		else:
			dpwm_score = 0
		if args.apwm: 
			apwm_score = im.score_pwm(aseq, re_apwm) * args.wapwm
		else:
			apwm_score = 0
		iscores[intron] = ilen_score + imm_score + dpwm_score + apwm_score
		dscores[intron] = dpwm_score
		ascores[intron] = apwm_score
	for exon in iso['exons']:
		iso['score'] += escores[exon]
	for intron in iso['introns']:
		iso['score'] += iscores[intron]
	iso['score'] -= len(iso['introns']) * args.icost * 100

abc_isoforms = sorted(abc_isoforms, key=lambda iso: iso['score'], reverse=True)
abc_isoforms = abc_isoforms[:args.limit]

'''
for a in abc_isoforms:
	if a['score'] != 0:
		print(a['beg'], a['end'], a['exons'], a['introns'], a['score'])
'''

iso_weights = []
iso_total = 0
for iso in abc_isoforms:
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
for iso in abc_isoforms:
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
print('# wb id:', wbgene)
print('# coordinates:', coor)
print('# strand:', strand)
print('# length:', len(seq))
print('# donors:', len(dons))
print('# acceptors:', len(accs))
print('# icost:', args.icost)
print('# trials:', trials)
print('# isoforms:', len(abc_isoforms))
print('# complexity:', f'{im.get_entropy(iso_probs):.4f}')

gff_writer = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
gff_writer.writerow([name, 'abc_isogen', 'gene', iso['beg']+1, iso['end']+1,
	'.', '+', '.', 'ID=Gene-' + name])
gff_writer.writerow([])
count = 0
for iso in abc_isoforms:
	if count <= args.limit - 1:
		iso_prob_f = '{:.5e}'.format(iso_probs[count])
		gff_writer.writerow([name, 'abc_isogen', 'mRNA', iso['beg']+1, 
			iso['end']+1, iso_prob_f, '+', '.', 'ID=iso-'+name+'-'+
			str(count+1)+';Parent=Gene-'+name])
		for exon in iso['exons']:
			escore_f = '{:.5e}'.format(escores[exon])
			efreq_f = '{:.5e}'.format(exon_freqs[exon])
			gff_writer.writerow([name, 'abc_isogen', 'exon', exon[0]+1,
				exon[1]+1, iso_prob_f, '+', '.', 'Parent='+'iso-'+name+'-'
				+str(count+1)+';score='+str(escore_f)+';exfreq='+str(efreq_f)])
		for intron in iso['introns']:
			iscore_f = '{:.5e}'.format(iscores[intron])
			ifreq_f = '{:.5e}'.format(intron_freqs[intron])
			gtscore_f = '{:5e}'.format(dscores[intron])
			agscore_f = '{:5e}'.format(ascores[intron])
			gff_writer.writerow([name, 'abc_isogen', 'intron', intron[0]+1,
				intron[1]+1, iso_prob_f, '+', '.', 'Parent='+'iso-'+name+'-'
				+str(count+1)+';score='+str(iscore_f)+';infreq='+str(ifreq_f)
				+';dscore='+str(gtscore_f)+';ascore='+str(agscore_f)])
		gff_writer.writerow([])
		count += 1



