import isosort_lib as isl
import argparse
import modelib as ml
# use ln -s for modelib in gff_analysis/

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
	help='fasta file')
parser.add_argument('wb_gff', type=str, metavar='<file>',
	help='gff file')
'''
parser.add_argument('--exon_len', type=str, metavar='<file>',
	help='exon length model .tsv')
parser.add_argument('--intron_len', type=str, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--exon_mm', type=str, metavar='<file>',
	help='exon markov model .tsv')
parser.add_argument('--intron_mm', type=str, metavar='<file>',
	help='intron markov model .tsv')
parser.add_argument('--donor_pwm', type=str, metavar='<file>',
	help='donor pwm .tsv')
parser.add_argument('--acceptor_pwm', type=str, metavar='<file>',
	help='donor pwm .tsv')
parser.add_argument('--icost', required=False, type=float, default=0.0, 
	metavar='<float>', help='intron cost %(default).2d')
'''
args = parser.parse_args()

seq = isl.get_seq(args.fasta)

wbginfo = isl.get_wbgene_info(args.wb_gff, seq)


print(wbginfo)






# need to include icost in wb gene score

args = parser.parse_args()


