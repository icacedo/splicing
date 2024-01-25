import isosort_lib as isl
import argparse
import os
import json

'''
parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
	help='wormbase fasta file for one gene')
parser.add_argument('wb_gff', type=str, metavar='<file>',
	help='wormbase nnotation gff file')
parser.add_argument('apcgen_gff', type=str, metavar='<file>',
	help='apc generated gff file')

args = parser.parse_args()

apcgen_isos = isl.amass_info(args.fasta, args.wb_gff, args.apcgen_gff)

for i in apcgen_isos:
	print(i, apcgen_isos[i])
'''

parser = argparse.ArgumentParser()
parser.add_argument('wb_dir', type=str, metavar='<directory>',
	help='directory with wormbase fasta/gff files')
parser.add_argument('apcgen_dir', type=str, metavar='<directory>',
	help='directory with apc generated gff files')

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
	isoforms_info = isl.amass_info(fa, wb_gff, apcgen_gff)
	jstr = json.dumps(isoforms_info, indent = 4)
	with open(f'out/{gID}.json', 'w') as outfile:
		outfile.write(jstr)








