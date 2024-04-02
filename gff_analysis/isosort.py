import isosort_lib as isl
import argparse
import os
import json

parser = argparse.ArgumentParser(
	description='adds labels to apc generated isoforms,' 
				' writes to .json files in out/')
parser.add_argument('wb_dir', type=str, metavar='<directory>',
	help='directory with wormbase fasta/gff files')
parser.add_argument('apcgen_dir', type=str, metavar='<directory>',
	help='directory with apc generated gff files')

parser.add_argument('--elen', type=str, metavar='<file>',
	help='exon length model .tsv')
parser.add_argument('--ilen', type=str, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--emm', type=str, metavar='<file>',
	help='exon length model .tsv')
parser.add_argument('--imm', type=str, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--dpwm', type=str, metavar='<file>',
	help='donor site pwm .tsv')
parser.add_argument('--apwm', type=str, metavar='<file>',
	help='acceptor site pwm .tsv')

parser.add_argument('--icost', required=False, type=float, default=22.0, 
	metavar='<float>', help='intron cost %(default).2d') 

args = parser.parse_args()

fastas = []
gffs = []
for fname in os.listdir(args.wb_dir):
	if fname.endswith('fa'): fastas.append(fname)
	if fname.endswith('gff3'): gffs.append(fname)

paths = {}
for fname in fastas:
	gID = fname.split('.')[1]
	paths[gID] = [args.wb_dir+fname]

for fname in gffs:
	gID = fname.split('.')[1]
	paths[gID] += [args.wb_dir+fname]

for fname in os.listdir(args.apcgen_dir):
	gID = fname.split('.')[1]
	paths[gID] += [args.apcgen_dir+fname]

os.makedirs('out/', exist_ok=True)

for gID in paths:
	fa = paths[gID][0]
	wb_gff = paths[gID][1]
	apcgen_gff = paths[gID][2]
	isoforms_info = isl.amass_info(fa, wb_gff, apcgen_gff, args.elen, 
									args.ilen, args.emm, args.imm, 
									args.dpwm, args.apwm, args.icost)
	jstr = json.dumps(isoforms_info, indent = 4)
	with open(f'out/{gID}.json', 'w') as outfile:
		outfile.write(jstr)








