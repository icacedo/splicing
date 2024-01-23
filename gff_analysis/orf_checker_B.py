# second version of orf_checker_lib.py
# reformat output dictionary
	
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
		wbg['mRNA'] = []
		wbg['exons'] = []
		for line in fp.readlines():
			line = line.rstrip()
			sline = line.split('\t')
			if sline[2] == 'mRNA':
				name = sline[0]+'-wb'
				wbg['mRNA'] = [int(sline[3]), int(sline[4])]
			if sline[2] == 'CDS':
				wbg['exons'].append((int(sline[3]), int(sline[4])))
		wbginfo[name] = wbg
		
	for gID in wbginfo:
		for ft in wbginfo[gID]:
			if ft == 'exons':
				wbginfo[gID][ft] = sorted(wbginfo[gID][ft]) 

	return wbginfo

def check_CDS(wbginfo):

	ntsum = 0
	for gID in wbginfo:
		for ft in wbginfo[gID]:
			if ft == 'exons':
				ntsum = 0
				for ex in wbginfo[gID][ft]:
					ntsum += ex[1] - ex[0] + 1
	if ntsum%3 == 0:
		wbginfo[gID]['in_frame'] = True
	else:
		wbginfo[gID]['in_frame'] = False

	return wbginfo
				
seq = get_seq(args.fasta)

wbginfo = get_wbgene_info(args.wb_gff, seq)

wbginfo = check_CDS(wbginfo)

print(wbginfo)

def get_apcgen_info(apcgen_gff):

	apc_isos = {}
	with open(apcgen_gff, 'r') as fp:
		icount = 0
		gID = ''
		for line in fp.readlines():
			line = line.rstrip()
			if 'name' in line:
				gID = line.split(' ')[2]
			sline = line.split('\t')
			if len(sline) < 9: continue
			if sline[2] == 'gene': continue
			if sline[2] == 'mRNA':
				icount += 1
				apc_isos[f'{gID}-{icount}'] = [sline]
			if sline[2] == 'exon':
				apc_isos[f'{gID}-{icount}'] += [sline]
			if sline[2] == 'intron':
				apc_isos[f'{gID}-{icount}'] += [sline]
	
	for iso in apc_isos:
		print(iso)
		mRNA = []
		exons = []
		introns = []
		for ft in apc_isos[iso]:
			print(ft)
			if ft[2] == 'mRNA':
				mRNA.append(ft[3])
				mRNA.append(ft[4])
			if ft[2] == 'exon':
				exons.append((ft[3], ft[4]))
			if ft[2] == 'intron':
				introns.append((ft[3], ft[4]))

		print(mRNA)	
		print(exons)
		print(introns)
		break
		
get_apcgen_info(args.apcgen_gff)
#apc = get_apcgen_info(args.apcgen_gff)
'''
for i in apc:
	for j in apc[i]:
		print(i, j)
	break
'''	










