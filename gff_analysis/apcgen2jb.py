import argparse

parser = argparse.ArgumentParser()
parser.add_argument('gff')

args = parser.parse_args()

# ch.9940 is on the negative strand, top iso matches wb
# WBGene0004483
# ch.9727 is on the positive strand, top iso matches wb
# WBGene00044472

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

print('# isoforms')

with open(args.gff, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		if len(line) < 9: continue
		if line[2] == 'gene': continue
			# remove gene group so isos can be selected individually
			#beg = coor[0] + int(line[3])
			#end = coor[0] + int(line[4])
			#gfs = (
			#	f'{chrom}\t{source}\tgene\t{beg}\t{end}\t.\t'
			#	f'{strand}\t.\t{line[8]}'
			#)
			#print(gfs)
		if line[2] == 'mRNA':
			line[0] = chrom	
			line[3] = str(coor[0] + int(line[3]))
			line[4] = str(coor[0] + int(line[4]))
			line[8] = line[8].split(';')[0]
			print('\t'.join(line))
		if line[2] == 'exon':
			if strand == '+':
				beg = coor[0] + int(line[3])
				end = coor[0] + int(line[4])
			if strand == '-':
				beg = coor[1] - int(line[4]) + 2
				end = coor[1] - int(line[3]) + 2
			par = line[8].split(';')[0]
			escore = escores[(int(line[3]), int(line[4]))]
			efs = (
				f'{chrom}\t{source}\texon\t{beg}\t{end}\t.\t'
				f'{strand}\t.\t{par};escore={escore}'
			)
			print(efs)

		if line[2] == 'intron':
			if strand == '+':
				beg = coor[0] + int(line[3])
				end = coor[0] + int(line[4])
			if strand == '-':
				beg = coor[1] - int(line[4]) + 2
				end = coor[1] - int(line[3]) + 2
			par = line[8].split(';')[0]
			iscore = iscores[(int(line[3]), int(line[4]))]
			dsite = intron[0]
			asite = intron[1]
			dscore = dscores[dsite]
			ascore = ascores[asite]
			ifs = ( 
				f'{chrom}\t{source}\tintron\t{beg}\t{end}\t.\t'
				f'{strand}\t.\t{par};iscore={iscore};dscore='
				f'{dscore};ascore={ascore}'
			)
			print(ifs)

print('# introns')

for intron in iscores:
	if strand == '+':
		beg = coor[0] + intron[0]
		end = coor[0] + intron[1]
	if strand == '-':
		beg = coor[1] - intron[1] + 2
		end = coor[1] - intron[0] + 2
	score = iscores[intron]
	ifs = (
		f'{chrom}\t{source}\tintron\t{beg}\t{end}\t{score}\t'
		f'{strand}\t.\tID={score}'
	)
	print(ifs)

print('# exons')

for exon in escores:
	if strand == '+':
		beg = coor[0] + exon[0]
		end = coor[0] + exon[1]
	if strand == '-':
		beg = coor[1] - exon[1] + 2
		end = coor[1] - exon[0] + 2
	score = escores[exon]
	efs = (
		f'{chrom}\t{source}\texon\t{beg}\t{end}\t{score}\t'
		f'{strand}\t.\tID={score}'
	)
	print(efs)

print('# donors')

for dsite in dscores:
	if strand == '+':
		beg = coor[0] + dsite
		end = beg + 4
	if strand == '-':
		end = coor[1] - dsite + 2
		beg = end - 4
	score = dscores[dsite]
	dfs = (
		f'{chrom}\t{source}\tdonor\t{beg}\t{end}\t{score}\t'
		f'{strand}\t.\tID={score}'
	)
	print(dfs)

print('# acceptors')

for asite in ascores:
	if strand == '+':
		end = coor[0] + asite
		beg = end - 5
	if strand == '-':
		beg = coor[1] - asite + 2
		end = beg + 5
	score = ascores[asite]
	afs = (
		f'{chrom}\t{source}\tacceptor\t{beg}\t{end}\t{score}\t'
		f'{strand}\t.\tID={score}'
	)
	print(afs)




	

