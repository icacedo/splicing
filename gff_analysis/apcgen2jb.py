import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gff')
parser.add_argument('fasta')

args = parser.parse_args()

escores = {}
iscores = {}
dscores = {}
ascores = {}
name = None
wbid = None
source = None
chrom = None
coor = None
strand = None
with open(args.gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if '# name' in line[0]:
			name = line[0].split(' ')[2]
		if '# wb id' in line[0]:
			wbid = line[0].split(' ')[3]
		if '# coordinates' in line[0]:
			loc = line[0].split(' ')[2]
			chrom = loc.split(':')[0]
			coor = loc.split(':')[1]
			coor = (int(coor.split('-')[0]), int(coor.split('-')[1]))
		if '# strand' in line[0]:
			strand = line[0].split(' ')[2]
		if len(line) == 9:
			if line[2] == 'gene':
				source = line[1]
			if line[2] == 'exon': 
				exon = (int(line[3]), int(line[4]))
				atbs = line[8].split(';')
				escore = float(atbs[1].split('=')[1])
				if exon not in escores:
					escores[exon] = escore
			if line[2] == 'intron':
				intron = (int(line[3]), int(line[4]))
				atbs = line[8].split(';')
				iscore = float(atbs[1].split('=')[1])
				dscore = float(atbs[3].split('=')[1])
				ascore = float(atbs[4].split('=')[1])
				if intron not in iscores:
					iscores[intron] = iscore
					dscores[intron[0]] = dscore
					ascores[intron[1]] = ascore

seq = ''
with open(args.fasta, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'): continue
		seq += line

for intron in iscores:
	beg = coor[0] + intron[0]
	end = coor[0] + intron[1]
	score = iscores[intron]
	print(f'{chrom}\t{source}\tintron\t{beg}\t{end}\t',
			f'{strand}\t.\tParent=Gene-{name};score={score}:name=jeans')

