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
			score = 0 if score == '.' else float(score)
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

	dist, details = isoform.expdiff(i1, i2)
	print(dist)
	for exon, p1, p2 in details:
		print(f'{exon[0]}\t{exon[1]}\t{p1:.6f}\t{p2:.6f}')




