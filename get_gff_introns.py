import sys

gff = sys.argv[1]

introns = []
with open(gff) as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if len(line) == 1: continue
		if line[2] == 'intron':
			iline = line[2:6]
			introns.append(' '.join(iline))
for i in introns:
	print(i)
print('*****')

introns = set(introns)

for i in introns:
	print(i)
print('*****')

info = []
with open(gff) as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if len(line) == 1: continue
		if line[1] == 'WormBase':
			wline = line[1:6]
			info.append(' '.join(wline))

# look at wormbase to find canonical isoform?
# can see what intron is a part of the wormbase entry
for wb in info:
	print(wb)




