import argparse
import isomod as im
import apc_model_lib as aml
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

seqid, seq = im.read_fasta(args.fasta)
seq_info = seqid.split(' ')
coor = seq_info[1]
strand = seq_info[2]
wbgene = seq_info[3].split(':')[1]
# 13301 is the shortest gene, quick to test
if args.gff:
	dons, accs = im.read_gff_sites(seq, args.gff) 
else:
	dons, accs = im.get_gtag(seq)

abc_isoforms, trials = im.abc(dons, accs, args.maxs, args.minin, 
							  args.minex, args.flank, seq)
'''
if args.elen:
	re_elen = im.read_len(args.elen)
else:
	re_elen = None
'''
re_elen = im.read_len(args.elen) if args.elen else None
re_ilen = im.read_len(args.elen) if args.ilen else None
re_emm = im.read_mm(args.emm) if args.emm else None
re_imm = im.read(args.imm) if args.imm else None
re_dpwm = im.read_pwm(args.dpwm) if args.dpwm else None
re_apwm = im.read_pwm(args.apwm) if args.apwm else None

# CONVERT PROBS TO SCORES
for iso in abc_isoforms:
	for exon in iso['exons']:
		elen_score = im.score_len(re_elen, exon)
		emm_score = im.score_mm(re_emm, exon, seq)
		print(emm_score)
		# finished getting score
	break

'''
# is the same
import apc_model_lib as aml
eseqs = ['ACTGATGCATGCATGC', 'GCTACGTA', 'GTCGCGTGTGACCCGAT']
mmsc, mmpb, order = aml.make_mm(eseqs)
print(dict(sorted(mmpb.items())))
print(dict(sorted(mmsc.items())))

import isoform as isf

mm = isf.create_markov(eseqs, 3, 0, 0)
print(dict(sorted(mm.items())))

print(math.log2(1/0.25), '#$##')
print('#####')
probs = [0.5, 0.1, 0.3, 0.1]
total = 0
for p in probs:
	total += math.log2(p/0.25)
print(total)
'''




'''
	total_iso_score = 0
	total_escore = 0
	total_iscore = 0
	for exon in iso['exons']:	

		total_escore += im.score_len(re_elen, exon) * args.welen
		print(im.score_len(re_elen, exon) * args.welen, '@@@')
		print(exon)
		total_escore += im.score_mm(re_emm, exon, seq) * args.wemm
		print(total_escore, '#####')
		exon_scores[exon] = total_escore
	for intron in iso['introns']:
		total_iscore += im.score_len(re_ilen, intron) * args.wilen
		total_iscore += im.score_mm(re_imm, intron, seq) * args.wimm
		dseq, aseq = im.get_daseq(intron, seq)
		dscore = im.score_pwm(dseq, re_dpwm) * args.wdpwm
		iscore = im.score_pwm(aseq, re_apwm) * args.wapwm
		total_iscore += dscore + iscore
		intron_scores[intron] = total_iscore
		gtag_scores[intron] = (dscore, iscore)
	total_iso_score = total_escore + total_iscore
	total_iso_score -= len(iso['introns']) * args.icost * 100
	iso['score'] = total_iso_score
	print(iso)


abc_isoforms = sorted(abc_isoforms, key=lambda iso: iso['score'], 
					  reverse=True)

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
			escore_f = '{:.5e}'.format(exon_scores[exon])
			efreq_f = '{:.5e}'.format(exon_freqs[exon])
			gff_writer.writerow([name, 'abc_isogen', 'exon', exon[0]+1,
				exon[1]+1, iso_prob_f, '+', '.', 'Parent='+'iso-'+name+'-'
				+str(count+1)+';score='+str(escore_f)+';exfreq='+str(efreq_f)])
		for intron in iso['introns']:
			iscore_f = '{:.5e}'.format(intron_scores[intron])
			ifreq_f = '{:.5e}'.format(intron_freqs[intron])
			gtscore_f = '{:5e}'.format(gtag_scores[intron][0])
			agscore_f = '{:5e}'.format(gtag_scores[intron][1])
			gff_writer.writerow([name, 'abc_isogen', 'intron', intron[0]+1,
				intron[1]+1, iso_prob_f, '+', '.', 'Parent='+'iso-'+name+'-'
				+str(count+1)+';score='+str(iscore_f)+';infreq='+str(ifreq_f)
				+';dscore='+str(gtscore_f)+';ascore='+str(agscore_f)])
		gff_writer.writerow([])
		count += 1
'''


