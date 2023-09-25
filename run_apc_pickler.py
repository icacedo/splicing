import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('apc_dir', type=str, metavar='<directory>',
	help='input directory with all apc fasta and gff files')
parser.add_argument('--read_gff', action='store_true', 
	help='get donor/acceptor sites from gff')

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

args = parser.parse_args()

fastas = {}
gffs = {}
for gfile in os.listdir(args.apc_dir):
	gfile = gfile.split('.')
	if gfile[2] == 'fa':
		fastas[gfile[1]] = '.'.join(gfile)
	if gfile[2] == 'gff3':
		gffs[gfile[1]] = '.'.join(gfile)

program = 'apc_pickler.py'
apc_dir = args.apc_dir
outdir = 'apc_pickles/'
os.makedirs(os.path.dirname(outdir), exist_ok=True)
max_splice = args.max_splice
min_intron = args.min_intron
min_exon = args.min_exon
flank = args.flank

if args.read_gff:
	print('wowo')
	for fID in fastas:
		fpath = apc_dir + fastas[fID]
		gpath = apc_dir + gffs[fID]
		print(fpath)
		print(gpath)
		subprocess.run(f'python3 {program} {fpath} --gff {gpath}'
			f' --max_splice {max_splice} --min_intron {min_intron}'
			f' --min_exon {min_exon} --flank {flank}', shell=True)
		break
else:
	for fID in fastas:
		fpath = apc_dir + fastas[fID]	
		subprocess.run(f'python3 {program} {fpath}'
			f' --max_splice {max_splice} --min_intron {min_intron}'
			f' --min_exon {min_exon} --flank {flank}', shell=True)
		break



