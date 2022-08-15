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
parser.add_argument('-gff3', type=str)
parser.add_argument('-fasta', type=str)

arg = parser.parse_args()

# gff reader
# move to seqlib when done
# APC dataset should only have + strands, ignore the few - strands

def gff_reader(gff3):

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
	
CDS_lines = gff_reader(arg.gff3)

CDS_regions = get_CDS_regions(CDS_lines)

CDS_sequences = {}

for ID, seq in seqlib.read_fasta(arg.fasta):
	ID = ID.split()
	ch_id = ID[0]
	CDS_seq = ''
	for CDS in CDS_regions[ch_id]:
		start = int(CDS[0])
		end = int(CDS[1])
		CDS_seq += seq[start-1:end]

	print(CDS_seq)
	CDS_sequences[ch_id] = CDS_seq
	
print(CDS_sequences)
		

		

		
































