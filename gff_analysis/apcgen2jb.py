import argparse
from grimoire.genome import Reader
from grimoire.feature import Feature

parser = argparse.ArgumentParser()
parser.add_argument('fasta')
parser.add_argument('gff')

args = parser.parse_args()

genome = Reader(gff=args.gff, fasta=args.fasta)
dna = next(genome)

loc, STR, gene = dna.desc.split()
chrom, coords = loc.split(':')
BEG, END = coords.split('-')
BEG = int(BEG)
END = int(END)

exos = []
ints = []
accs = []
dons = []
for f in dna.ftable.features:
	coor = (f.beg, f.end)
	if f.type == 'intron':
		if coor not in ints: ints.append(coor)
		if f.beg not in dons: dons.append(f.beg)
		if f.end not in accs: accs.append(f.end)
	elif f.type == 'exon':
		if coor not in exos: exos.append(coor)

escores = {}
iscores = {}
dscores = {}
ascores = {}
with open(args.gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if len(line) == 9:
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

gffs = []	
for beg, end in exos:
	gffs.append(Feature(dna, beg, end, STR, 'exon', source='model',
				score=escores[(beg, end)]))
for beg, end in ints:
	s = iscores[(beg, end)]
	gffs.append(Feature(dna, beg, end, STR, 'intron', source='model',
				score=iscores[(beg, end)]))
for d in dons:
	gffs.append(Feature(dna, d, d+1, STR, 'donor', source='model',
				score=dscores[d]))
for a in accs:
	gffs.append(Feature(dna, a, a+1, STR, 'acceptor', source='model',
				score=ascores[a]))
for f in gffs:
	dna.ftable.add_feature(f)

if STR == '-': dna.revcomp()

for f in dna.ftable.features:
	col9 = ''
	if f.id: col9 = f'ID={f.id}'
	if f.pid:
		pids = ','.join(f.pid)
		if col9: col9 += f';Parent={pids}'
		else: col9 = f'Parent={pids}'

	print('\t'.join((chrom, f.source, f.type,
			str(f.beg + BEG), str(f.end + END),
			str(f.score), f.strand, '.', col9)))




