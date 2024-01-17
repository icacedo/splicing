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
				stCDS[f'nts-{count}'] = cdsnt
			# check if all CDS together is in frame	
			# should always be in frame for wormbase			
			if ntsum%3 == 0: stCDS['in frame?'] = 'yes'
			if ntsum%3 != 0: stCDS['in frame?'] = 'no'
			wbgS[gID] = stCDS
		
			return wbgS

def get_apciso_info(apcgen_gff):
	
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

	aisosinfos = {}
	for apc_iso in apc_isos:
		aisoinfo = {}
		ecount = 0
		icount = 0
		for info in apc_isos[apc_iso]:
			if info[2] == 'mRNA':
				name = info[0] + '-' + str(apc_iso)
				aisoinfo['mRNA'] = (int(info[3]), int(info[4]))
				aisoinfo['prob'] = float(info[5])
			if info[2] == 'exon':
				ecount += 1
				aisoinfo['exon'+f'-{ecount}'] = (int(info[3]), int(info[4]))
			if info[2] == 'intron':
				icount += 1
				aisoinfo['intron'+f'-{icount}'] = (int(info[3]), int(info[4]))
		aisosinfos[name] = aisoinfo

	return aisosinfos

seq = get_seq(args.fasta)

wbg_infos = get_wbgene_info(args.wb_gff, seq)

apcgen_isos = get_apciso_info(args.apcgen_gff)

for gn in wbg_infos:
	clist = []
	for ft in wbg_infos[gn]:
		if 'CDS' in ft:
			clist.append(ft)
	wbstart = wbg_infos[gn][clist[0]][0]
	wbstop = wbg_infos[gn][clist[-1]][1]
print(wbstart, wbstop)

print(wbg_infos)
for iso in apcgen_isos:
	print(iso, apcgen_isos[iso])
	break

for isoid in apcgen_isos:
	exons = []
	ntsum = 0
	for ft in apcgen_isos[isoid]:
		if ft == 'mRNA':
			mRNA = apcgen_isos[isoid][ft]
		if 'exon' in ft:
			exons.append(ft)
	print(isoid, mRNA, (wbstart, wbstop))
	first = apcgen_isos[isoid][exons[0]]	
	print(exons[0], seq[wbstart-1:wbstart+2], seq[first[1]-3:first[1]], wbstart, first[1])
	ntsum += first[1] - wbstart + 1
	for ex in exons:
		if ex == exons[0] or ex == exons[-1]: continue
		beg = apcgen_isos[isoid][ex][0]
		end = apcgen_isos[isoid][ex][1]
		print(ex, seq[beg-1:beg+2], seq[end-3:end], beg, end)
		ntsum += end - beg + 1
	last = apcgen_isos[isoid][exons[-1]]
	print(exons[-1], seq[last[0]-1:last[0]+2], seq[wbstop-3:wbstop], last[0], wbstop)
	ntsum += wbstop - last[0] + 1
	print(ntsum/3, 'is in frame?', wbstop, last[0])
	print('### NEXT ISO ###')

# need to test an example that has only one CDS





















