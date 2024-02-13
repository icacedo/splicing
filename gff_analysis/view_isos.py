import sys
import isosort_lib as isl
import json
import os
import re

fasta = sys.argv[1]
jfile = sys.argv[2]

seq = isl.get_seq(fasta)
seq = seq.lower()

seq_sites = re.sub('gt', 'GT', seq)
seq_sites = re.sub('ag', 'AG', seq_sites, flags=re.IGNORECASE)
seq_stops = re.sub('tag', 'TAG', seq_sites, flags=re.IGNORECASE)
seq_stops = re.sub('taa', 'TAA', seq_stops, flags=re.IGNORECASE)
seq_mod = re.sub('tga', 'TGA', seq_stops, flags=re.IGNORECASE)

def make_wb_sym(jfile, seq):

	name = os.path.basename(jfile).split('.')[0]
	with open(jfile, 'r') as jf: info = json.load(jf)
	
	endf = info[f'ch.{name}-wb']['exons'][0][0]
	fseq1 = seq[:endf-1]
	fsym1 = ''
	for i in range(len(fseq1)):
		fsym1 += '#'
	begf = info[f'ch.{name}-wb']['exons'][-1][1]
	fseq2 = seq[begf:]
	fsym2 = ''
	for i in range(len(fseq2)):
		fsym2 += '#'
	esyms = []
	for e in info[f'ch.{name}-wb']['exons']:
		eseq = seq[e[0]-1:e[1]]
		esym = ''
		for i in range(len(eseq)): esym += '*'
		esyms.append(esym)
	isyms = []
	for i in info[f'ch.{name}-wb']['introns']:
		iseq = seq[i[0]-1:i[1]]
		isym = ''
		for i in range(len(iseq)): isym += '-'
		isyms.append(isym)
		
	sym_seq_wb = ''
	sym_seq_wb += fsym1
	for i in range(len(esyms)):
		sym_seq_wb += esyms[i]
		if i > len(isyms)-1: continue
		sym_seq_wb += isyms[i]
	sym_seq_wb += fsym2

	return sym_seq_wb

sym_seq_wb = make_wb_sym(jfile, seq_mod)

for i in range(round(len(seq_mod)/80)):
	print(seq_mod[i*80:i*80+80])
	print(sym_seq_wb[i*80:i*80+80])
	print('')
