import sys

# ch.7157 ATG lines up with wormbase
# ch.9018 ATG lines up with wormbase
# five_prime_UTR in gff is where the start codon is 

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

with open(gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		sline = line.split('\t')
		#print(sline)
		if sline[2] == 'five_prime_UTR':
			print(sline)
			start = int(sline[4])
			print(seq[start:start+10], start) 

'''
for i in range(0, len(seq)):
	if seq[i:i+3] == 'ATG': 
		print(seq[i:i+3], i)

print(seq[126:136])
'''
