# goal is to swap the gene_id identfier position in 1pct.gffread.gtf
# to see if works with subread
# this did not make gffread converted file usable by subread

import sys
import re

f = sys.argv[1]
fp = open(f)

f2 = sys.argv[2]
fp2 = open(f2, 'w')

#count = 0
for line in fp.readlines():
	line = line.split('\t')
	string = line[8]
	strings = re.split('(;)', string)
	transcript_id = strings[0] + strings[1]
	try:
		gene_id = strings[2].strip() + strings[3]
	except IndexError: 
		gene_id = strings[2].strip() + ';'
	c9 = gene_id + ' ' + transcript_id
	fp2.write('\t'.join(line[0:8]) + '\t' + c9 + '\n')
	#count += 1
	#if count == 10: sys.exit()
	

