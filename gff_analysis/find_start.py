import sys

# use ch.7157

gff = sys.argv[1]
fasta = sys.argv[2]

with open(gff, 'r') as fp:
	for line in fp.readlines():
		print(line.rstrip())

with open(fasta, 'r') as fp:
	seq = ''
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			ID = line
		else:
			seq += line

for i in range(0, len(seq)):
	if seq[i:i+3] == 'ATG': 
		print(seq[i:i+3], i)

print(seq[355:598])

