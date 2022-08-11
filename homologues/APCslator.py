# goal is to find the proteins in the APC dataset
# make own translator-glue CDS's together from apc_test
# do a blastp protein to protein
# download all the protein sequences
# blast the APC proteomes together

import argparse
import seqlib

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>')
parser.add_argument('gff', type=str, metavar='<file>')

arg = parser.parse_args()
'''
for ID, seq in seqlib.read_fasta('input/apc_test/ch.10010.fa'):
	print(ID)
	print(seq)
'''	
# gff reader
# move to seqlib when done

fp = open(arg.gff)
for line in fp.readlines():
	line = line.rstrip()
	line = line.split()
	seqid = line[0]
	source = line[1]
	type = line[2]
	start = line[3]
	end = line[4]
	score = line[5]
	strand = line[6]
	phase = line[7]
	if source == 'WormBase':
		attributes = line[8]
	if type == 'CDS':
		print(line)
		print(attributes)
