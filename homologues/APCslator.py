# goal is to find the proteins in the APC dataset
# make own translator-glue CDS's together from apc_test
# do a blastp protein to protein
# download all the protein sequences
# blast the APC proteomes together

import argparse
import seqlib
import os
import sys

parser = argparse.ArgumentParser()
#parser.add_argument('--gff3', type=str)
#parser.add_argument('--fasta', type=str)
parser.add_argument('--path_to_apc', type=str)

arg = parser.parse_args()

# gff reader
# move to seqlib when done
# APC dataset should only have + strands, ignore the few - strands

def get_CDS_lines(gff3):

	fp = open(gff3)
	CDS_lines = []
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split()
		SEQID = line[0]
		SOURCE = line[1]
		TYPE = line[2]
		START = line[3]
		END = line[4]
		SCORE = line[5]
		STRAND = line[6]
		PHASE = line[7]
		if STRAND == '-': continue
		if SOURCE == 'WormBase':
			ATTRIBUTES = line[8]
		if TYPE == 'CDS':
			CDS_lines.append(line)
	
	return CDS_lines

def get_CDS_regions(CDS_lines):

	CDS_regions = {}
	starts_n_ends = []	
	for line in CDS_lines:
		startend = line[3], line[4]
		starts_n_ends.append(startend)
		CDS_regions[line[0]] = starts_n_ends

	return CDS_regions

def get_CDS_sequences(CDS_regions, fasta):

	CDS_sequences = {}
	for ID, seq in seqlib.read_fasta(fasta):
		ID = ID.split()
		ch_id = ID[0]
		CDS_seq = ''
		for CDS in CDS_regions[ch_id]:
			start = int(CDS[0])
			end = int(CDS[1])
			CDS_seq += seq[start-1:end]

		CDS_sequences[ch_id] = CDS_seq
	
	return CDS_sequences

# can put multiple sequences in the same file for command line blast
# also, do i need to worry about the phase?
# used this command for all gff files in the apc/ directory
# awk '//{print $8}' data/apc/*.gff3 | sort --unique
# there only '.'
# i will just run blastx for now, it doesn't need a specific start site


		


path = arg.path_to_apc

ext_gff = 'gff3'
ext_fasta = 'fa'
already_seen = []
for files in os.listdir(path):
	split_file = os.path.splitext(files)
	chID = split_file[0]
	if chID in already_seen: continue
	else: 
		already_seen.append(chID)
	gff_path = path + chID + '.' + ext_gff
	fasta_path = path + chID + '.' + ext_fasta  	
	CDS_lines = get_CDS_lines(gff_path)
	CDS_regions = get_CDS_regions(CDS_lines)
	print(get_CDS_sequences(CDS_regions,fasta_path))
	
	
'''	
filename = 'c_elegans.APC.CDS.fa'

for ID in CDS_sequences:
	seq = CDS_sequences[ID]
	with open(filename, 'w') as fo:	
		fo.write(ID+'.CDS'+'\n')
		lines = []
		for i in range(0, len(seq), 80):
			line = seq[i:i+80]
			lines.append(line)
		fo.write('\n'.join(lines))
'''































