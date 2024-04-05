import argparse
import isosort_lib as isl
import csv
import os
import json

parser = argparse.ArgumentParser()
parser.add_argument('fa_dir', type=str, metavar='<directory>',
	help='directory with wormbase fasta files')
parser.add_argument('j_dir', type=str, metavar='<directory>',
	help='directory with apc isoform json files')
parser.add_argument('g_list', type=str, metavar='<file>',
	help='csv file with list of genes with bad isos to examine')

args = parser.parse_args()

paths = {}
with open(args.g_list, 'r') as csvfile:
	reader = csv.DictReader(csvfile)
	for row in reader:
		ID = row['iso id'].split('.')[1]
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
	with open(jpath, 'r') as jf: info = json.load(jf)
	seq = isl.get_seq(fpath)
	seq_mod = isl.mod_seq(seq)	
	sym_seq_wb = isl.make_wb_sym(jpath, seq)
	frame_sym = isl.make_frame_sym(sym_seq_wb)
	with open(f'txt_views/{ID}.txt', 'w') as outfile:
		outfile.write(f'   {frame_sym}\n')
		outfile.write(f'   {seq_mod}\n')
		outfile.write(f'wb {sym_seq_wb}\n')
		wb_escores = info[f'ch.{ID}-wb']["escores"]
		wb_iscores = info[f'ch.{ID}-wb']["iscores"]
		wb_gtag_scores = info[f'ch.{ID}-wb']["gtag_scores"]
		es_line = isl.string_hyphs(seq, f'   exon scores: {wb_escores}')
		is_line = isl.string_hyphs(seq, f'   intron scores: {wb_iscores}')
		ga_line = isl.string_hyphs(seq, f'   gtag scores: {wb_gtag_scores}')
		outfile.write(f'{es_line}\n')
		outfile.write(f'{is_line}\n')
		outfile.write(f'{ga_line}\n')
		count = 1
		for iso in info:
			if iso == f'ch.{ID}-wb': continue
			sym_seq_apc = isl.make_apc_sym(iso, info, seq)
			wbin = info[f'ch.{ID}-wb']['introns']
			isoin = info[iso]['introns']
			if wbin == isoin:
				outfile.write(f'wb {sym_seq_apc}\n')
			else:
				if len(str(count)) == 1:
					outfile.write(f'{count}  {sym_seq_apc}\n')
				if len(str(count)) == 2:
					outfile.write(f'{count} {sym_seq_apc}\n')
			count += 1
			es_line = isl.string_hyphs(seq, 
								f'   exon scores: {info[iso]["escores"]}')
			is_line = isl.string_hyphs(seq, 
								f'   intron scores: {info[iso]["iscores"]}')
			ga_line = isl.string_hyphs(seq, 
								f'   gtag scores: {info[iso]["gtag_scores"]}')
			outfile.write(f'{es_line}\n')
			outfile.write(f'{is_line}\n')
			outfile.write(f'{ga_line}\n')
