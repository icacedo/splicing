import argparse
import os

parser = argparse.ArgumentParser(description=
	'writes file with command line arguments to be used in parallize')
parser.add_argument('apc_dir', type=str, metavar='<str>', 
	help='directory with apc fasta and gff files')
parser.add_argument('--read_gff', action='store_true', 
	help='get don/acc sites from gff files in apc dir')
parser.add_argument('--outfile', required=True, type=str, metavar='<str>', 
	help='name and directory of apc cmd file')
parser.add_argument('--gff_out', required=True, type=str, metavar='<str>',
	help='output directory for gff files')
parser.add_argument('--gff_name', required=True, type=str, metavar='<str>',
	help='name for gff files, e.g., gID.{name}.gff')

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
parser.add_argument('--icost', required=False, type=float, default=21,
	metavar='<float>', help='intron cost %(default).2d')
	
args = parser.parse_args()

apc_list = []
for fname in os.listdir(args.apc_dir):
	apc_list.append(args.apc_dir + fname)

fa_gff_pairs = {}
for fpath in sorted(apc_list):
	ID = fpath.split('.')[-2]
	if ID not in fa_gff_pairs:
		fa_gff_pairs[ID] = [fpath]
	else:
		fa_gff_pairs[ID] += [fpath]

os.makedirs(args.gff_out, exist_ok=True)

f = open(args.outfile, 'w')
for gID in fa_gff_pairs:
	fa_path = fa_gff_pairs[gID][0]
	gff_path = fa_gff_pairs[gID][1]
	gff_in = ''
	if args.read_gff:
		gff_in = ' --gff ' + gff_path
	f.write(
		f'python3 apc_isogen.py {fa_path}'
		f'{gff_in}'
		f' --min_intron {args.min_intron} --min_exon {args.min_exon}'
		f' --flank {args.flank} --limit {args.limit}'
		f' --exon_len {args.exon_len} --intron_len {args.intron_len}'
		f' --exon_mm {args.exon_mm} --intron_mm {args.intron_mm}'
		f' --donor_pwm {args.donor_pwm} --acceptor_pwm {args.acceptor_pwm}'
		f' --icost {args.icost}'
		f' > {args.gff_out}ch.{gID}.{args.gff_name}.gff'
		f'\n'
	)
f.close()






