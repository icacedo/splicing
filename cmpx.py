import argparse
import isoform


def get_introns(gff):
	introns = {}
	total = 0
	with open(gff) as fp:
		for line in fp.readlines():
			if line.startswith('#'): continue
			f = line.split()
			if len(f) < 8: continue
			if f[2] != 'intron': continue
			if f[6] != '+': continue
			beg, end, score = f[3], f[4], f[5]
			beg = int(beg)
			end = int(end)
			score = float(score)
			if (beg,end) not in introns: introns[(beg,end)] = 0
			introns[(beg,end)] += score
			total += score

	# convert to histogram
	for k in introns: introns[k] /= total

	return introns

if __name__ == '__main__':

	## Command Line Interface ##

	parser = argparse.ArgumentParser(
		description='Compare Expression')
	parser.add_argument('gff1', type=str, metavar='<file>',
		help='input gff file 1')
	parser.add_argument('gff2', type=str, metavar='<file>',
		help='input gff file 2')
	arg = parser.parse_args()

	i1 = get_introns(arg.gff1)
	i2 = get_introns(arg.gff2)

	# ensure all introns are in both collections
	for k in i1:
		if k not in i2: i2[k] = 0
	for k in i2:
		if k not in i1: i1[k] = 0

	p1 = []
	p2 = []
	for k in i1:
		p1.append(i1[k])
		p2.append(i2[k])

	print(isoform.manhattan(p1, p2))


