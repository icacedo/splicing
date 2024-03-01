import argparse
import isosort_lib as isl
import os

parser = argparse.ArgumentParser()
parser.add_argument('fa_dir', type=str, metavar='<directory>',
	help='directory with wormbase fasta files')
parser.add_argument('j_dir', type=str, metavar='<directory>',
	help='directory with apc isoform json files')
parser.add_argument('g_list', type=str, metavar='<file>',
	help='text file with list of bad isos to examine')

args = parser.parse_args()

paths = {}

with open(args.g_list, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('#'): continue
		ID = line.split('.')[1]
		paths[ID] = []

for fname in os.listdir(args.fa_dir):
	if fname.endswith('.gff3'): continue
	ID = fname.split('.')[1]
	if ID not in paths: continue	
	paths[ID].append(args.fa_dir+fname)

for fname in os.listdir(args.j_dir):
	ID = fname.split('.')[0]
	if ID not in paths: continue
	paths[ID].append(args.j_dir+fname)

os.makedirs('txt_views/', exist_ok=True)

for ID in paths:
	fpath = paths[ID][0]
	jpath = paths[ID][1]
	seq = isl.get_seq(fpath)
	seq_mod = isl.mod_seq(seq)
	sym_seq_wb = isl.make_wb_sym(jpath, seq)
	sym_seq_apc = isl.make_apc_sym(jpath, seq)
	frame_sym = isl.make_frame_sym(sym_seq_wb)
	isl.add_scores(jpath)
	break
	with open(f'txt_views/{ID}.txt', 'w') as outfile:
		outfile.write(f'{frame_sym}\n')
		outfile.write(f'{sym_seq_wb}\n')
		outfile.write(f'{seq_mod}\n')
		outfile.write(f'{sym_seq_apc}\n')

