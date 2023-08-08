# log of intron sequences/ total sequences
# use gff to get intron seqs
# calculate one intron cost for every isoform

import sys
import modelib as ml
import os

seq_alt = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'

gff = sys.argv[1]
fasta = sys.argv[2]
apc_dir = sys.argv[3]

def total_iseqs(gff, fasta):

	seq = None
	for seqid, seq in ml.read_fastas(fasta):
		seq = seq

	with open(gff, 'r') as fp:
		total_iseqs = 0
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split()
			if line[2] == 'intron':
				beg = int(line[3]) - 1 
				end = int(line[4])
				iseq = seq[beg:end]
				total_iseqs += len(iseq)
		return len(seq), total_iseqs

seq_len, iseqs_len = total_iseqs(gff, fasta)

print(seq_len, iseqs_len)

print(apc_dir)

gffs = []
fastas = []
for filename in os.listdir(apc_dir):
	f_path = os.path.join(apc_dir, filename)
	if f_path.endswith('.gff3'):
		gffs.append(f_path)
	if f_path.endswith('.fa'):
		fastas.append(f_path)

fa_gff_pairs = {}
for gff_path in gffs:
	gid = gff_path.split('.')[1]
	for fa_path in fastas:
		fid = fa_path.split('.')[1]
		if gid == fid:
			fa_gff_pairs[gff_path] = fa_path

for gff_key in fa_gff_pairs:
	total_bp, intron_bp = total_iseqs(gff_key, fa_gff_pairs[gff_key])
	print(total_bp, intron_bp)
	

