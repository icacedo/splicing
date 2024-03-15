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

with open(args.gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if len(line) < 9: continue
		if line[2] == 'gene':
			beg = coor[0] + int(line[3])
			end = coor[1] + int(line[4])
			print(f'{chrom}\t{source}\tgene\t{beg}\t{end}\t.\t',
					f'{strand}\t.\t{line[8]}')
		if line[2] == 'mRNA':
			line[0] = chrom	
			line[3] = str(coor[0] + int(line[3]))
			line[4] = str(coor[1] + int(line[4]))
			print('\t'.join(line))
		if line[2] == 'exon':
			beg = coor[0] + int(line[3])
			end = coor[0] + int(line[4])
			par = line[8].split(';')[0]
			print(f'{chrom}\t{source}\texon\t{beg}\t{end}\t.\t',
					f'{strand}\t.\t{par}')
		if line[2] == 'intron':	
			beg = coor[0] + int(line[3])
			end = coor[0] + int(line[4])
			par = line[8].split(';')[0]
			print(f'{chrom}\t{source}\tintron\t{beg}\t{end}\t.\t',
					f'{strand}\t.\t{par}')

#print('#')

for intron in iscores:
	beg = coor[0] + intron[0]
	end = coor[0] + intron[1]
	score = iscores[intron]
	print(f'{chrom}\t{source}\tintron\t{beg}\t{end}\t{score}\t',
			f'{strand}\t.\tID={score}')

#print('#')

for exon in escores:
	beg = coor[0] + exon[0]
	end = coor[0] + exon[1]
	score = escores[exon]
	print(f'{chrom}\t{source}\texon\t{beg}\t{end}\t{score}\t',
			f'{strand}\t.\tID={score}')

#print('#')

for dsite in dscores:
	beg = coor[0] + dsite
	end = beg + 4
	score = dscores[dsite]
	print(f'{chrom}\t{source}\tdonor\t{beg}\t{end}\t{score}\t',
			f'{strand}\t.\tID={score}')

#print('#')

for asite in ascores:
	end = coor[0] + asite
	beg = end - 5
	score = ascores[asite]
	print(f'{chrom}\t{source}\tacceptor\t{beg}\t{end}\t{score}\t',
			f'{strand}\t.\tID={score}')





	

