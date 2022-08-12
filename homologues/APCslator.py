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
#parser.add_argument('fasta', type=str)
#parser.add_argument('gff3', type=str)
parser.add_argument('path', type=str)

arg = parser.parse_args()
'''
for ID, seq in seqlib.read_fasta('input/apc_test/ch.10010.fa'):
	print(ID)
	print(seq)
'''	
# gff reader
# move to seqlib when done
# not sure what what make of + or - strand
# would i need to make the complement of the sequence in the fasta file?
# to get the CDS sequence on the - strand?
'''
fp = open(arg.gff3)
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
	if SOURCE == 'WormBase':
		ATTRIBUTES = line[8]
	if TYPE == 'CDS':
		print(line)
		print(ATTRIBUTES)
'''		

# *just ignore the - strands*
# are there any CDS regions on the - strand? YES
# in the apc dataset
# 1896 CDS regions on the + strand, and 14 on the - strand
# are there - strand CDS with + strand CDS in same sequence/geneID? YES

#countplus = 0
#countnega = 0
for item in os.listdir(arg.path):
	if item.endswith('gff3'):
		fp = open(arg.path+item)
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split()
			if line[2] == 'CDS':
				if line[6] == '+':
					print(line[0], line[2], line[3], line[4], line[6])
					#print(line[2],line[6])
					#countplus += 1
				if line[6] == '-':
					print(line)
					#print(line[2],line[3],line[4],line[6])
					#countnega += 1
#print(countplus, countnega)

# save output of this code to > strands
# use to find - strands and see all other strands:
# grep --color=always -e '^' -e '-' strands
































