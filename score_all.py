import argparse
import modelib as ml

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
# probabilistic models
parser.add_argument('--exon_len', required=False, type=int, metavar='<file>', 
	help='exon length model .tsv')
parser.add_argument('--intron_len', required=False, type=int, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--exon_mm', required=False, type=int, metavar='<file>',
	help='exon markov model .tsv')
parser.add_argument('--intron_mm', required=False, type=int, metavar='<file>',
	help='intron markov model .tsv')
parser.add_argument('--donor_pwm', required=False, type=int, metavar='<file>',
	help='donor pwm .tsv')
parser.add_argument('--acceptor_pwm', required=False, type=int, metavar='<file>',
	help='acceptor pwm .tsv')

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

apc_isoforms, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)

#exon_len_pdf, exon_len_sco = ml.read_exin_len(args.exon_len)
#intron_len_pdf, intron_len_sco = ml.read_exin_len(args.intron_len)
# need option for gff don/acceptor sites
#for iso in apc_isoforms:
#	d_seqs, a_seqs = get_donacc_seqs
#	print(iso)

# test geniso in arch/
# python3 geniso --min_intron 3 --min_exon 4 --flank 5 ../test_seq.fa
# https://useast.ensembl.org/info/website/upload/gff.html
