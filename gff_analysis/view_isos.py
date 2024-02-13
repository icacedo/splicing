import sys
import isosort_lib as isl
import json
import os
import re

fasta = sys.argv[1]
jfile = sys.argv[2]

seq = isl.get_seq(fasta)

seq_sites = re.sub('GT', 'gt', seq)
seq_sites = re.sub('AG', 'ag', seq_sites, flags=re.IGNORECASE)

print(seq_sites)

name = os.path.basename(jfile).split('.')[0]
with open(jfile, 'r') as jf:
	info = json.load(jf)	
	print(info[f'ch.{name}-wb'])
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
	
sym_seq = ''
sym_seq += fsym1
for i in range(len(esyms)):
	sym_seq += esyms[i]
	if i > len(isyms)-1: continue
	sym_seq += isyms[i]
sym_seq += fsym2

print(sym_seq)

print(len(seq_sites))
print(len(sym_seq))

print(1080/80)
print(round(1080/80))
for i in range(round(len(seq_sites)/80)):
	#print(len(seq_sites[i*80:i*80+80]))
	print(seq_sites[i*80:i*80+80])
	print(sym_seq[i*80:i*80+80])
	print('')
