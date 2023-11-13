import sys

gff = sys.argv[1]


isos = {
	'iso1': {
		'exons': [(101, 285), (355, 597)],
		'introns': [(286, 354)]
	}
}

print(isos)

with open(gff, 'r') as fp:
	count = 0
	for line in fp.readlines():
		line = line.rstrip()
		sline = line.split('\t')
		if len(sline) < 9: continue
		if sline[2] == 'mRNA':
			count += 1
			print('start', count)
		
		
			




		


