import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gff_1', type=str, metavar='<file>')
parser.add_argument('gff_2', type=str, metavar='<file>')

args = parser.parse_args()

def get_gff_intron_probs(gff):

	with open(gff) as fp:
		introns = {}
		total_score = 0
		for line in fp.readlines():
			if line.startswith('#'): continue
			line = line.rstrip()
			line = line.split('\t')
			if len(line) < 8: continue
			if line[2] == 'intron' and line[6] == '+':
				beg = int(line[3])
				end = int(line[4])
				score = line[5]
				intron = (beg, end)
				if score == '.': introns[intron] = 0
				else: 
					if intron not in introns: introns[intron] = 0
					introns[intron] += float(score)
					total_score += float(score)
		
		for i in introns:
			introns[i] = introns[i]/total_score

		return introns
			
introns1 = get_gff_intron_probs(args.gff_1)
introns2 = get_gff_intron_probs(args.gff_2)

for i in introns1:
	if i not in introns2:
		introns2[i] = 0
for i in introns2:
	if i not in introns1:
		introns1[i] = 0

introns1 = dict(sorted(introns1.items(), 
				key=lambda item: item[1], reverse=True))

dd = 0
for i in introns1:
	print(i, introns1[i], introns2[i])
	d = introns1[i] - introns2[i]
	print(abs(d))
	dd += abs(d)

print(dd)
