# goal is to swap the gene_id identfier position in 1pct.gffread.gtf
# to see if works with subread

import sys
import re

f = sys.argv[1]
fp = open(f)

count = 0
for line in fp.readlines():
	line = line.split('\t')
	string = line[8]
	#print(string)
	strings = re.split('(;)', string)
	print(strings)
	#transcript_id = strings[0] + strings[1]
	#gene_id = strings[2] + strings[3]
	#print(transcript_id)
	#print(gene_id) 
	count += 1
	if count == 2: sys.exit()
