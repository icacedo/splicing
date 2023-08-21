import sys

gff = sys.argv[1]

introns = []
with open(gff) as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if len(line) == 1: continue
		if line[2] == 'intron':
			iline = line[2:5]
			introns.append(' '.join(iline))

introns = set(introns)

for i in introns:
	print(i)




