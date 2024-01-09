import sys

# five_prime_UTR in gff is where the start codon is 
# ch.5202 does not have a five_prime_UTR feature 
# use ch.216 to test
# INTRON DOES NOT NEED TO BE IN FRAME!


gff = sys.argv[1]
fasta = sys.argv[2]

with open(fasta, 'r') as fp:
	seq = ''
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			ID = line
		else:
			seq += line

# use CDS instead of 5' and 3' UTR
# 
with open(gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		sline = line.split('\t')
		# need to handle an exeption when there is no 5' utr
		# or just use CDS
		if sline[2] == 'CDS':
			print(sline)
		if sline[2] == 'five_prime_UTR':
			print(sline)
			# for 216, 5' utr ends at 111, gene starts at 112
			start = int(sline[4])
			print(seq[start:start+3], seq[start], start, '***')
		if sline[2] == 'three_prime_UTR':
			print(sline)
			end = int(sline[3]) - 2
			print(seq[end-2:end+1], seq[end], end, '@@@')
			print(start, end)
			print((end-start-1+2)/3)

s = 'ATGCGGATAG'
print(s[0], s[9])

# everything is in frame
for i in range(start, end, 3):
	print(seq[i:i+3])
	
