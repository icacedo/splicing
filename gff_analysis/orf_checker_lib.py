import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
	help='wormbase fasta file for one gene')
parser.add_argument('wb_gff', type=str, metavar='<file>',
	help='wormbase nnotation gff file')
parser.add_argument('apcgen_gff', type=str, metavar='<file>',
	help='apc generated gff file')

args = parser.parse_args()

def get_seq(fasta):

	with open(fasta, 'r') as fp:
		seq = ''
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('>'):
				gID = line
			else:
				seq += line

	return seq	
	
def get_wbgene_info(wb_gff, seq):
	
	wbginfo = {}
	with open(wb_gff, 'r') as fp:
		wbg = {}
		count = 0
		for line in fp.readlines():
			line = line.rstrip()
			sline = line.split('\t')
			if sline[2] == 'mRNA':
				name = sline[0]+'-wb'
				wbg['mRNA'] = (int(sline[3]), int(sline[4]))
			if sline[2] == 'CDS':
				count += 1
				wbg['CDS'+f'-{count}'] = (int(sline[3]), int(sline[4]))
		wbginfo[name] = wbg
		
		wbgS = {}
		for gID in wbginfo:
			cdsnames = []
			for ft in wbginfo[gID]:
				if ft == 'mRNA':
					mRNA = wbginfo[gID][ft]
				if 'CDS' in ft:
					cdsnames.append(ft)
			cdscrs = []
			for cn in cdsnames:
				cdscrs.append(wbginfo[gID][cn])
			stCDS = {}
			stCDS['mRNA'] = mRNA
			count = 0
			ntsum = 0
			for csc in sorted(cdscrs):
				count += 1 
				stCDS['CDS'+f'-{count}'] = csc
				bcds = seq[csc[0]-1:csc[0]+2]
				ecds = seq[csc[1]-3:csc[1]]
				cdsnt = csc[1] - csc[0] + 1
				ntsum += cdsnt
				# first and last codons in the CDS
				stCDS[f'codons-{count}'] = (bcds, ecds)
				# number of nucleotides in the CDS
				#stCDS[f'nts-{count}'] = cdsnt
			# check if all CDS together is in frame	
			# should always be in frame for wormbase			
			if ntsum%3 == 0: stCDS['in frame?'] = 'yes'
			if ntsum%3 != 0: stCDS['in frame?'] = 'no'
			wbgS[gID] = stCDS
		
			return wbgS

# creates dictionary for each isoform in the gff
def get_apcgen_iso_info(apcgen_gff, wbstart, wbstop):
	
	apc_isos = {}
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
				apc_isos[count] = [sline]
			else:
				apc_isos[count] += [sline]

	agen_isos = {}
	for apc_iso in apc_isos:
		agen_iso = {}
		ecount = 0
		icount = 0
		exons = []
		ntsum = 0
		for info in apc_isos[apc_iso]:
			if info[2] == 'mRNA':
				name = info[0] + '-' + str(apc_iso)
				agen_iso['mRNA'] = (int(info[3]), int(info[4]))
				agen_iso['prob'] = float(info[5])
				agen_iso['wbstartstop'] = (wbstart, wbstop)
			if info[2] == 'exon':
				ecount += 1
				agen_iso[f'exon-{ecount}'] = (int(info[3]), int(info[4]))
				exons.append((int(info[3]), int(info[4])))
			if info[2] == 'intron':
				icount += 1
				agen_iso[f'intron-{icount}'] = (int(info[3]), int(info[4]))
		ecount2 = 1
		first = exons[0]
		last = exons[-1]
		c0 = seq[wbstart-1:wbstart+2]
		c1 = seq[first[1]-3:first[1]]
		agen_iso[f'codons-{ecount2}'] = (c0, c1)
		ntsum += first[1] - wbstart + 1
		for ex in exons:
			if ex == exons[0] or ex == exons[-1]: continue	
			ecount2 += 1
			c1 = seq[ex[0]-1:ex[0]+2]
			c2 = seq[ex[1]-3:ex[1]]
			agen_iso[f'codons-{ecount2}'] = (c1, c2)
			ntsum += ex[1] - ex[0] + 1
		last = exons[-1]
		ecount2 += 1
		c3 = seq[last[0]-1:last[0]+2]
		c4 = seq[wbstop-3:wbstop]
		agen_iso[f'codons-{ecount2}'] = (c3, c4) 
		ntsum += wbstop - last[0] + 1
		agen_isos[name] = agen_iso
		if ntsum%3 == 0:
			agen_iso['in_frame'] = 'yes'
		if ntsum%3 != 0:
			agen_iso['in_frame'] = 'no'
	
	return agen_isos

# checks if the number of exons is the same as wormbase
# adds info to apcgen_isos
# exons for apcgen, CDS for wb
def check_exon_count(apcgen_isos, wbg_info):

	count = 0
	for iso in apcgen_isos:
		gid = iso.split('-')[0] + '-wb'
		ccount = 1
		cnum = 0
		for wft in wbg_info[gid]:
			if wft == f'CDS-{ccount}':
				ccount += 1
				cnum += 1
		ecount = 1
		enum = 0
		for aft in apcgen_isos[iso]:	
			if aft == f'exon-{ecount}':
				ecount += 1
				enum += 1
		if enum != cnum:
			apcgen_isos[iso]['dif_ex_num'] = 'yes'
		if enum == cnum:
			apcgen_isos[iso]['dif_ex_num'] = 'no'
	
	return apcgen_isos

# get start and stop codon positions from wb
def get_wb_start_stop(wbg_info):
	for gn in wbg_info:
		clist = []
		for ft in wbg_info[gn]:
			if 'CDS' in ft:
				clist.append(ft)
		wbstart = wbg_info[gn][clist[0]][0]
		wbstop = wbg_info[gn][clist[-1]][1]
	
	return wbstart, wbstop



seq = get_seq(args.fasta)

wbg_info = get_wbgene_info(args.wb_gff, seq)

wbstart, wbstop = get_wb_start_stop(wbg_info)

apcgen_isos = get_apcgen_iso_info(args.apcgen_gff, wbstart, wbstop)

print(wbg_info)
print('#####')
apcgen_isos = check_exon_count(apcgen_isos, wbg_info)

# check if same frame as wb
count = 0
for iso in apcgen_isos:
	if apcgen_isos[iso]['dif_ex_num'] == 'yes':
		apcgen_isos[iso]['wb_frame'] = 'no'
		print(apcgen_isos[iso])
	if apcgen_isos[iso]['dif_ex_num'] == 'no':
		gid = iso.split('-')[0] + '-wb'
		wb_exlist = []
		ccount = 1
		for wft in wbg_info[gid]:
			if wft == f'CDS-{ccount}':
				wb_exlist.append(wbg_info[gid][wft])
				ccount += 1
		#for aft in apcgen_isos[iso]:
		#	if 'exon
		print(wb_exlist)
				
		print(gid)
		
			
	if count == 2: break
	count += 1



print('#####')
see = 0
for i in apcgen_isos:
	print(apcgen_isos[i])
	if see == 2: break
	see += 1


# need to test an example that has only one CDS
# will always have at least 2 CDS, bc need at least one intron
# but will APC generate isos with no introns?
# ch.11934 has 3 CDS
# ch.216 has 2 CDS
# ch.4738 has 4 CDS 
# ch.4741, 2nd isoform is given an extra exon/intron




















