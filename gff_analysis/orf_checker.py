import sys

fasta = sys.argv[1]
wb_gff = sys.argv[2]
apcgen_gff = sys.argv[3]

with open(fasta, 'r') as fp:
	seq = ''
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			gID = line
		else:
			seq += line

# find wormbase start/stop and orf
# ch.11934 has 3 CDS
# CDS appear both earlier to later and later to earlier
# how can i reorder the CDS?

wbg = {}
with open(wb_gff, 'r') as fp:
	wbginfo = {}
	count = 0
	for line in fp.readlines():
		line = line.rstrip()
		sline = line.split('\t')	
		if sline[2] == 'mRNA':
			print(sline)
			name = sline[0]+'-wb'
			wbginfo['mRNA'] = (int(sline[3]), int(sline[4]))
		if sline[2] == 'CDS':
			print(sline)
			count += 1
			wbginfo['CDS'+f'-{count}'] = (int(sline[3]), int(sline[4]))
	wbg[name] = wbginfo
print('#####')

print('*WBG*', wbg)
wbgS = {}
for k1 in wbg:
	cdsnames = []
	for k2 in wbg[k1]:
		if k2 == 'mRNA':
			mRNA = wbg[k1][k2]
		if 'CDS' in k2:
			cdsnames.append(k2)
	cdscrs = []
	for cn in cdsnames:
		cdscrs.append(wbg[k1][cn])
	# rewrite dict w sorted CDS
	stCDS = {}
	print(wbgS, '******')
	#print('name: ', k1)
	#print('mRNA: ', mRNA)
	stCDS['mRNA'] = mRNA
	count = 0
	for csc in sorted(cdscrs):
		count += 1 
		stCDS['CDS'+f'-{count}'] = csc
		#print('CDS'+f'-{count}: ', csc)
	wbgS[k1] = stCDS
print('*WBG SORTED*', wbgS)

for gn in wbg:
	clist = []
	for ft in wbg[gn]:
		if 'CDS' in ft:
			clist.append(ft)
	start = wbg[gn][clist[0]][0]
	end = wbg[gn][clist[-1]][1]
	print('wb start and stop codons', start, end)


# goal: print ATG, # of bases, STOP	
# stop codons: TAG, TAA, TGA

aisos = {}
with open(apcgen_gff, 'r') as fp:
	count = 0
	for line in fp.readlines():
		line = line.rstrip()
		sline = line.split('\t')
		if len(sline) < 9: continue
		if sline[2] == 'gene': 
			continue
		if sline[2] == 'mRNA':
			count += 1
			aisos[count] = [sline]
		else:
			aisos[count] += [sline]

isosinfos = {}
for aiso in aisos:
	isoinfo = {}
	ecount = 0
	icount = 0
	for info in aisos[aiso]:
		if info[2] == 'mRNA':
			name = info[0] + '-' + str(aiso)
			isoinfo['mRNA'] = (int(info[3]), int(info[4]))
			isoinfo['prob'] = float(info[5])
		if info[2] == 'exon':
			ecount += 1
			isoinfo['exon'+f'-{ecount}'] = (int(info[3]), int(info[4]))
		if info[2] == 'intron':
			icount += 1
			isoinfo['intron'+f'-{icount}'] = (int(info[3]), int(info[4]))
	isosinfos[name] = isoinfo
	# testing with only the top isoform
	break
print(isosinfos)

print('##########')

# now compare sequences

# apc isoforms need to use wb ATG
for inm in isosinfos:
	print(inm)

'''
			print(gn, ft, wbg[gn][ft])
			start = wbg[gn][ft][0] - 1
			print(start, seq[start:start+3])
			stop = wbg[gn][ft][1]
			print(stop, seq[stop-3:stop])

print('***')
seek = 'ATGNGCGCGCGNTAG'
print(len(seek))
print(seek[0:3])
print(seek[12:15])
print((len(seek)-6)/3)
print(seek[3:12])
'''
# need to test an example that has only one CDS

