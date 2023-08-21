import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gff_1', type=str, metavar='<file>')
parser.add_argument('gff_2', type=str, metavar='<file>')

args = parser.parse_args()

def get_gff_introns(gff):

	with open(gff) as fp:

		ilines = []
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			if len(line) < 9: continue
			if line[2] == 'intron':
				iline = line[2:5]
				ilines.append(' '.join(iline))

		return ilines

introns1 = get_gff_introns(args.gff_1)
introns2 = get_gff_introns(args.gff_2)
print(ilines)


