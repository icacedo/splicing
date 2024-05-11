import argparse
import os

parser = argparse.ArgumentParser(description=
	'writes file with command line arguments to be used in parallize')
parser.add_argument('apc_dir', type=str, metavar='<str>', 
	help='directory with apc fasta and gff files')
parser.add_argument('--weights', required=False, type=str, metavar='<str>',
	help='file with individual weights for each gene')
parser.add_argument('--read_gff', action='store_true', 
	help='get don/acc sites from gff files in apc dir')
parser.add_argument('--outfile', required=True, type=str, metavar='<str>', 
	help='name and directory of apc cmd file')
parser.add_argument('--gff_out', required=True, type=str, metavar='<str>',
	help='output directory for gff files')
parser.add_argument('--gff_name', required=True, type=str, metavar='<str>',
	help='name for gff files, e.g., gID.{name}.gff')

# apc parameters
parser.add_argument('--maxs', required=False, type=int, default=3,
	metavar='<int>', help='maximum number of splicing events %(default)d')
parser.add_argument('--minin', required=False, type=int, default=35,
	metavar='<int>', help='minimum length of intron %(default)d')
parser.add_argument('--minex', required=False, type=int, default=25,
	metavar='<int>', help='minimum length of exon %(default)d')
parser.add_argument('--flank', required=False, type=int, default=100,
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

if args.weights:
	ftwts = {}
	with open(args.weights, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			fit = line[0]
			wdpwm = line[1]
			wapwm = line[2]
			wemm = line[3]
			wimm = line[4]
			welen = line[5]
			wilen = line[6]
			icost = line[7]
			gid = line[8]
			ftwts[gid] = {
				'fit': fit,
				'wdpwm': wdpwm,
				'wapwm': wapwm,
				'wemm': wemm,
				'wimm': wimm,
				'welen': welen,
				'wilen': wilen,
				'icost': icost
			}

elen = f' --elen {args.elen}' if args.elen else ''
ilen = f' --ilen {args.ilen}' if args.ilen else ''
emm = f' --emm {args.emm}' if args.emm else ''
imm = f' --imm {args.imm}' if args.imm else ''
dpwm = f' --dpwm {args.dpwm}' if args.dpwm else ''
apwm = f' --apwm {args.apwm}' if args.apwm else ''

os.makedirs(args.gff_out, exist_ok=True)

f = open(args.outfile, 'w')
for gID in fa_gff_pairs:
	gid = f'ch.{gID}'
	fa_path = fa_gff_pairs[gID][0]
	gff_path = fa_gff_pairs[gID][1]
	gff_in = ''
	if args.read_gff:
		gff_in = ' --gff ' + gff_path
	wts = ''
	if args.weights:
		wts = (
			f' --wdpwm {ftwts[gid]['wdpwm']} --wapwm {ftwts[gid]['wapwm']} '
			f'--wemm {ftwts[gid]['wemm']} --wimm {ftwts[gid]['wimm']} '
			f'--welen {ftwts[gid]['welen']} --wilen {ftwts[gid]['wilen']} '
			f'--icost {ftwts[gid]['icost']}'
		)
	f.write(
		f'python3 apc_isogen.py {fa_path}'
		f'{gff_in}'
		f'{wts}'
		f' --maxs {args.maxs} --minin {args.minin} --minex {args.minex}'
		f' --flank {args.flank} --limit {args.limit}'
		f'{elen}' 
		f'{ilen}'
		f'{emm}'
		f'{imm}'
		f'{dpwm}'
		f'{apwm}'
		f' > {args.gff_out}ch.{gID}.{args.gff_name}.gff'
		f'\n'
	)
f.close()






