import argparse
import modelib as ml

# params for test seq: maxs 100, minin 3, minex 4, flank 5
parser = argparse.ArgumentParser(
	description='generate and score alternative isoforms')
parser.add_argument('fasta', type=str, metavar='<file>',
	help='input single sequence fasta file')
parser.add_argument('--gff', type=str, metavar='<file>', required=False,
	help='input .gff3 for single gene')

# apc parameters
parser.add_argument('--max_splice', required=False, type=int, default=100,
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
if args.intron_len:
	re_ilen_pdf, re_ilen_log2 = ml.read_exin_len(args.intron_len)
if args.exon_mm:
	re_emm_prob, re_emm_log2 = ml.read_exin_mm(args.exon_mm)
if args.intron_mm:
	re_imm_prob, re_imm_log2 = ml.read_exin_mm(args.intron_mm)
if args.donor_pwm:
	re_dppm, re_dpwm = ml.read_pwm(args.donor_pwm)
if args.acceptor_pwm:
	re_appm, re_apwm = ml.read_pwm(args.acceptor_pwm)

# there are differences in how mm and pwm are calculated in geniso
# need to investigate further
for iso in apc_isoforms:
	#print(iso)
	exon_lengths, intron_lengths = ml.get_exin_lengths(iso)
	elen_score = ml.get_len_score(exon_lengths, re_elen_log2)
	ilen_score = ml.get_len_score(intron_lengths, re_ilen_log2)
	print('len ex in', elen_score, ilen_score)	
	exon_seqs, intron_seqs = ml.get_exin_seqs(iso, seq)
	emm_score = ml.get_mm_score(exon_seqs, re_emm_log2)
	imm_score = ml.get_mm_score(intron_seqs, re_imm_log2)
	print('mm ex in', emm_score, imm_score)
	donor_seqs, acceptor_seqs = ml.get_donacc_seqs(iso, seq)
	dpwm_score = ml.get_pwm_score(donor_seqs, re_dpwm)
	apwm_score = ml.get_pwm_score(acceptor_seqs, re_apwm)
	print('pwm d a', dpwm_score, apwm_score)
	iso['score'] = elen_score + ilen_score + emm_score + imm_score + \
		dpwm_score + apwm_score
	print(iso)
	break

# log of intron seq/total seq
# this is for all the sequences
# there is the same cost for each isoform
# test geniso in arch/
# python3 geniso --min_intron 3 --min_exon 4 --flank 5 ../test_seq.fa
# https://useast.ensembl.org/info/website/upload/gff.html
print('##########')
# test old apc code
import isoform_fixed as isof

name, seq = next(isof.read_fasta(args.fasta))
txs, info = isof.all_possible(seq, args.min_intron, args.min_exon,
	args.max_splice, args.flank)

for tx in txs:
	print(tx)
print('##########')
for iso in apc_isoforms:
	print(iso)

# testing geniso
'''
python3 geniso_test sqtst.fa --min_intron 3 --min_exon 4 --flank 5 --dpwm arch/data/donor.pwm --apwm arch/data/acceptor.pwm --emm arch/data/exon.mm --imm arch/data/intron.mm --elen arch/data/exon.len --ilen arch/data/intron.len 
'''

