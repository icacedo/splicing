import argparse
import modelib as ml
import csv
import sys

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

exon_len_score = None
intron_len_score = None
exon_mm_score = None
intron_mm_score = None
donor_pwm = None
acceptor_pwm = None

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

for iso in apc_isoforms:
	exon_lengths, intron_lengths = ml.get_exin_lengths(iso)
	elen_score = ml.get_len_score(exon_lengths, re_elen_log2, ea, eb, eg)
	ilen_score = ml.get_len_score(intron_lengths, re_ilen_log2, ia, ib, ig)
	#print('len ex in', elen_score, ilen_score)	
	exon_seqs, intron_seqs = ml.get_exin_seqs(iso, seq)
	emm_score = ml.get_mm_score(exon_seqs, re_emm_log2)
	imm_score = ml.get_mm_score(intron_seqs, re_imm_log2, 'GT', 'AG')
	#print('mm ex in', emm_score, imm_score)
	donor_seqs, acceptor_seqs = ml.get_donacc_seqs(iso, seq)
	dpwm_score = ml.get_pwm_score(donor_seqs, re_dpwm)
	apwm_score = ml.get_pwm_score(acceptor_seqs, re_apwm)
	#print('pwm d a', dpwm_score, apwm_score)
	score = elen_score + ilen_score + emm_score + imm_score + \
		dpwm_score + apwm_score
	#print(iso['introns'], len(iso['introns']))
	score -= len(iso['introns']) * args.icost
	iso['score'] = score
	#print(iso)

# print as gff
# https://useast.ensembl.org/info/website/upload/gff.html

seqname = seqid
source = 'apc_isogen'

apc_isoforms = sorted(apc_isoforms, key=lambda iso: iso['score'], reverse=True)

weights = []
total = 0
for iso in apc_isoforms:
	weight = 2 ** iso['score']
	weights.append(weight)
	total += weight
	
probs = []
for w in weights:
	probs.append(w / total)

# still differences in geniso calculated probabilities
# differences by a few decimals

# to calc an icost:
# when choosing an intron cost, just start from an arbitrary number like 0-100
# and see how the cost affects the outcome

print('# seqid: ' + seqid)
print('# length: ' + str(len(seq)))
print('# donors: ' + str(len(dons)))
print('# acceptors: ' + str(len(accs)))
print('# trials: ' + str(trials))
print('# isoforms: ' + str(len(apc_isoforms)))
print('# complexity: ' + 'later')

name = seqid.split(' ')[0]

# is everything on the + strand? and not worried about frame?
gff_writer = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
gff_writer.writerow([name, 'apc_isogen', 'gene', iso['beg']+1, iso['end']+1,
	 '.', '+', '.', 'ID=Gene-' + name])
gff_writer.writerow([])
count = 0
for iso in apc_isoforms:
	if count <= args.limit - 1:
		probs_f = '{:.5e}'.format(probs[count])
		gff_writer.writerow([name, 'apc', 'mRNA', iso['beg']+1, 
			iso['end']+1, probs_f, '+', '.', 'ID=iso-'+name+'-'+str(count+1)+
				';Parent=Gene-'+name])
		for exon in iso['exons']:
			gff_writer.writerow([name, 'apc', 'exon', exon[0]+1, exon[1]+1, 
				probs_f, '+', '.', 'Parent='+'iso-'+name+'-'+str(count+1)])
		for intron in iso['introns']:
			gff_writer.writerow([name, 'apc', 'intron', intron[0]+1, intron[1]+1,
				probs_f, '+', '.', 'Parent='+'iso-'+name+'-'+str(count+1)])
		gff_writer.writerow([])
		count += 1

	





# test geniso in arch
# python3 geniso --min_intron 3 --min_exon 4 --flank 5 ../test_seq.fa
# https://useast.ensembl.org/info/website/upload/gff.html

'''
For the test seq 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
and first isoform
parameters for testing: min_intron 3, min_exon 4, flank 5
{'seq': 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA', 'beg': 5, 'end': 42, 
'exons': [(5, 8), (17, 42)], 'introns': [(9, 16)], 'score': 0}
individual scores should be (there are small decimal differences between this
output and geniso) using reformatted length models from geniso:

score_all:
len ex in -102.85725982788392 -100.0
mm ex in -2.9090535606612224 -0.7277110496653949
pwm d a -0.8923380000000001 -7.053903999999999
{'seq': 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA', 'beg': 5, 'end': 42, 
'exons': [(5, 8), (17, 42)], 'introns': [(9, 16)], 'score': -214.44026643821053}

geniso:
apwm -7.053965734900418
dpwm -0.8923498962049816
elen -102.85725982788392
ilen -100.0
emm -2.9090604877439743
imm -0.7014542894184398
{'seq': 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA', 'beg': 5, 'end': 42, 
'exons': [(5, 8), (17, 42)], 'introns': [(9, 16)], 'score': -214.41409023615174}
'''

