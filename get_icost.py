# log of intron sequences/ total sequences
# use gff to get intron seqs
# calculate one intron cost for every isoform

# problem: the number of intron sequnces is ~2x the total number of bp
# can't use every reported intron

# solution? limit counted intron seqs to only WormBase introns

import sys
import modelib as ml
import os
import math

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

def total_iseqs_canon(gff, fasta):

	seq = None
	for seqid, seq in ml.read_fastas(fasta):
		seq = seq

	with open(gff, 'r') as fp:
		total_iseq_len = 0
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split()
			if line[1] == 'WormBase' and line[2] == 'intron':
				beg = int(line[3]) - 1
				end = int(line[4])
				iseq_len = end - beg
				total_iseq_len += iseq_len
		return len(seq), total_iseq_len

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

# some genes have more than one canonical intron
# they are also counted

intron = 0
total = 0
for gff_key in fa_gff_pairs:
	total_bp, intron_bp = total_iseqs_canon(gff_key, fa_gff_pairs[gff_key])
	if gff_key.split('/')[3] == 'ch.61.gff3': print('##################################')
	print(total_bp, intron_bp, gff_key.split('/')[3])
	intron += intron_bp
	total += total_bp

print(intron, total)

icost = math.log2(intron/total)
print(icost)
print('#####')

with open(gff, 'r') as fp:
	total_iseq_len = 0
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split()
		if line[1] == 'WormBase' and line[2] == 'intron':
			print(line[1:6])
			beg = int(line[3]) - 1
			end = int(line[4])
			iseq_len = end - beg
			print(iseq_len)

# i don't remember what this was for
'''
print('**************************')

saved_scores = []
with open(gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split()
		if line[2] == 'intron':
			if line[5] == '.':
				saved_scores.append(0)
			else:
				saved_scores.append(float(line[5]))
print(max(saved_scores))
'''
