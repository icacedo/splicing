import re
import os
import json
import apc_model_lib as aml

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
		wbg['introns'] = []
		for line in fp.readlines():
			line = line.rstrip()
			sline = line.split('\t')
			if sline[2] == 'mRNA':
				name = sline[0]+'-wb'
				wbg['mRNA'] = [int(sline[3]), int(sline[4])]
				WBGene = sline[8].split(':')[2]
				wbg['Parent=Gene'] = WBGene
			if sline[2] == 'CDS':
				wbg['exons'].append((int(sline[3]), int(sline[4])))
			if sline[1] == 'WormBase' and sline[2] == 'intron':
				wbg['introns'].append((int(sline[3]), int(sline[4])))
		wbginfo[name] = wbg
		
	for gID in wbginfo:
		for ft in wbginfo[gID]:
			if ft == 'exons':
				wbginfo[gID][ft] = sorted(wbginfo[gID][ft]) 

	return wbginfo

# read in .tsv files for probabilistic models
def score_wb_iso(seq, wbginfo, elen, ilen, emm, imm, dpwm, apwm, icost):
		
	re_elen_pdf, re_elen_log2 = aml.read_exin_len(elen)
	ea, eb, eg = aml.read_len_params(elen)

	re_ilen_pdf, re_ilen_log2 = aml.read_exin_len(ilen)
	ia, ib, ig = aml.read_len_params(ilen)

	re_emm_prob, re_emm_log2 = aml.read_exin_mm(emm)
	re_imm_prob, re_imm_log2 = aml.read_exin_mm(imm)

	re_dppm, re_dpwm = aml.read_pwm(dpwm)
	re_appm, re_apwm = aml.read_pwm(apwm)

	for gene in wbginfo:
		wbginfo[gene]['escores'] = []
		wbginfo[gene]['total_icost_score'] = 0
		for exon in wbginfo[gene]['exons']:
			# include region before start site and after 100 bp flank
			if exon == wbginfo[gene]['exons'][0]: 
				exon = (101, exon[1])
			if exon == wbginfo[gene]['exons'][-1]:
				exon = (exon[0], len(seq)-100)
			exon = (exon[0]-1, exon[1]-1) # adjust indexing
			elen_score = aml.get_exin_len_score(exon, re_elen_log2, ea, eb, eg)
			emm_score = aml.get_exin_mm_score(exon, seq, re_emm_log2)
			escore = elen_score + emm_score
			escore = float('{:.5e}'.format(escore))
			wbginfo[gene]['total_icost_score'] += escore
			wbginfo[gene]['escores'].append(escore)
		wbginfo[gene]['iscores'] = []
		wbginfo[gene]['gtag_scores'] = []
		for intron in wbginfo[gene]['introns']:
			intron = (intron[0]-1, intron[1]-1) # adjust indexing
			ilen_score = aml.get_exin_len_score(intron, re_ilen_log2, ia, ib, ig)
			imm_score = aml.get_exin_mm_score(intron, seq, re_imm_log2, 'GT', 'AG')
			dseq, aseq = aml.get_donacc_seq(intron, seq)
			dpwm_score = aml.get_donacc_pwm_score(dseq, re_dpwm)
			apwm_score = aml.get_donacc_pwm_score(aseq, re_apwm)
			wbginfo[gene]['gtag_scores'].append((dpwm_score, apwm_score))
			iscore = ilen_score + imm_score + dpwm_score + apwm_score
			iscore = float('{:.5e}'.format(iscore))
			wbginfo[gene]['iscores'].append(iscore)
			wbginfo[gene]['total_icost_score'] += iscore
		wbginfo[gene]['total_icost_score'] -= len(wbginfo[gene]['introns']) \
																	* icost
		return wbginfo

def check_CDS(info):

	ntsum = 0
	for gID in info:	
		ntsum = 0
		for ex in info[gID]['exons']:
			ntsum += ex[1] - ex[0] + 1
		if ntsum%3 == 0:
			info[gID]['in_frame'] = True
		else:
			info[gID]['in_frame'] = False
	
	return info

def get_start_stop(wbginfo):
	
	for gn in wbginfo:
		clist = []
		wbstart = wbginfo[gn]['exons'][0][0]	
		wbstop = wbginfo[gn]['exons'][-1][1]

	return wbstart, wbstop	
	
def get_apcgen_info(seq, apcgen_gff, wbstart, wbstop, dpwm, apwm):

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
	
	re_dppm, re_dpwm = aml.read_pwm(dpwm)
	re_appm, re_apwm = aml.read_pwm(apwm)

	apcgen_isos = {}	
	for iso in apc_isos:
		apcgen_isos[iso] = {}
		mRNA = []
		prob = float
		escores = []
		iscores = []
		gtag_scores = []
		exons = []
		introns = []
		for ft in apc_isos[iso]:
			if ft[2] == 'mRNA':
				mRNA.append(ft[3])
				mRNA.append(ft[4])
				prob = float(ft[5])
			if ft[2] == 'exon':
				exons.append((int(ft[3]), int(ft[4])))
				escore = ft[8].split(';')[1]
				escore = float(escore.split('=')[1])
				escores.append(escore)
			if ft[2] == 'intron':
				dsite = int(ft[3]) - 1
				asite = int(ft[4]) - 1
				dseq, aseq = aml.get_donacc_seq((dsite, asite), seq)
				dpwm_score = aml.get_donacc_pwm_score(dseq, re_dpwm)
				apwm_score = aml.get_donacc_pwm_score(aseq, re_apwm)
				dpwm_score = float('{:.5e}'.format(dpwm_score))
				apwm_score = float('{:.5e}'.format(apwm_score))	
				gtag_scores.append((dpwm_score, apwm_score))
				introns.append((int(ft[3]), int(ft[4])))
				iscore = ft[8].split(';')[1]
				iscore = float(iscore.split('=')[1])
				iscores.append(iscore)
		exons2 = []
		first = (wbstart, exons[0][1])
		exons2.append(first)
		for ex in exons:
			if ex == exons[0]: continue
			if ex == exons[-1]: continue
			exons2.append(ex)
		last = (exons[-1][0], wbstop)
		exons2.append(last)
		apcgen_isos[iso]['mRNA'] = mRNA
		apcgen_isos[iso]['prob'] = prob
		apcgen_isos[iso]['exons'] = exons2
		apcgen_isos[iso]['escores'] = escores
		apcgen_isos[iso]['introns'] = introns
		apcgen_isos[iso]['iscores'] = iscores
		apcgen_isos[iso]['gtag_scores'] = gtag_scores

	return apcgen_isos

def check_exon_count(apcgen_isos, wbg_info):

	for iso in apcgen_isos:
		gID = iso.split('-')[0] + '-wb'
		wb_enum = len(wbg_info[gID]['exons'])
		apc_enum = len(apcgen_isos[iso]['exons'])
		apcgen_isos[iso]['dif_exon'] = apc_enum - wb_enum 

	return apcgen_isos

def check_wb_frame(apcgen_isos, wbg_info):

	for iso in apcgen_isos:
		if apcgen_isos[iso]['dif_exon'] != 0:
			apcgen_isos[iso]['wb_frame'] = False
		if apcgen_isos[iso]['dif_exon'] == 0:
			gID = iso.split('-')[0] + '-wb'
			wb_ex = wbg_info[gID]['exons']
			apc_ex = apcgen_isos[iso]['exons']
			if wb_ex == apc_ex: 
				apcgen_isos[iso]['wb_frame'] = True
			else:
				apcgen_isos[iso]['wb_frame'] = False
	
	return apcgen_isos

def get_codons(apcgen_isos, seq):

	for iso in apcgen_isos:
		codons = []
		for ex in apcgen_isos[iso]['exons']:
			first = seq[ex[0]-1:ex[0]+2]
			last = seq[ex[1]-3:ex[1]]
			codons.append([first, last])
		apcgen_isos[iso]['codons'] = codons
	
	return apcgen_isos

def find_PTCs(apcgen_isos, seq):

	for iso in apcgen_isos:
		if apcgen_isos[iso]['wb_frame'] == False:
			CDS = ''
			for ex in apcgen_isos[iso]['exons']:
				eseq = seq[ex[0]-1:ex[1]]
				CDS += eseq
			shift = 0
			PTCs = []
			for i in range(len(CDS)):
				codon = CDS[i+shift:i+shift+3]
				if len(codon) == 3: 
					if codon in ['TAG', 'TAA', 'TGA']:
						PTCs.append((i+shift+1, codon))	
				shift += 2
			if len(PTCs) > 0:
				apcgen_isos[iso]['PTC'] = PTCs
			else:
				apcgen_isos[iso]['PTC'] = False
		else:
			apcgen_isos[iso]['PTC'] = False
		
	return apcgen_isos	

def amass_info(fasta, wb_gff, apcgen_gff, elen, 
				ilen, emm, imm, dpwm, apwm, icost):
	
	seq = get_seq(fasta)
	wbg_info = get_wbgene_info(wb_gff, seq)
	wbg_info = check_CDS(wbg_info)
	wbg_info = get_codons(wbg_info, seq)
	wbg_info = score_wb_iso(seq, wbg_info, elen, ilen, 
							emm, imm, dpwm, apwm, icost)
	wbstart, wbstop = get_start_stop(wbg_info)
	apcgen_isos = get_apcgen_info(seq, apcgen_gff, wbstart, wbstop, dpwm, apwm)
	apcgen_isos = check_CDS(apcgen_isos)
	apcgen_isos = check_exon_count(apcgen_isos, wbg_info)
	apcgen_isos = check_wb_frame(apcgen_isos, wbg_info)
	apcgen_isos = get_codons(apcgen_isos, seq)
	apcgen_isos = find_PTCs(apcgen_isos, seq)
	apcgen_isos.update(wbg_info)

	return apcgen_isos

# ch.4738 has a short first exon 
# all wb genes have at least 1 intron
# but will APC generate isos with no introns?
# ch.11934 has 3 CDS, ordered last to first
# ch.216 has 2 CDS
# ch.4738 has 4 CDS 
# ch.4741, 2nd isoform is given an extra exon/intron
# ch.241 has 3 CDS, 2nd isoform has correct firt and last exons
# but middle exon is cut out as an intron, then goes out of frame

# section to write txt iso view files

def mod_seq(seq):
	
	seq = seq.lower()
	seq_sites = re.sub('gt', 'GT', seq)
	seq_sites = re.sub('ag', 'AG', seq_sites, flags=re.IGNORECASE)
	seq_stops = re.sub('tag', 'TAG', seq_sites, flags=re.IGNORECASE)
	seq_stops = re.sub('taa', 'TAA', seq_stops, flags=re.IGNORECASE)
	seq_mod = re.sub('tga', 'TGA', seq_stops, flags=re.IGNORECASE)

	return seq_mod	

def make_wb_sym(jfile, seq):

	name = os.path.basename(jfile).split('.')[0]
	with open(jfile, 'r') as jf: info = json.load(jf)
	
	endf = info[f'ch.{name}-wb']['exons'][0][0]
	fseq1 = seq[:endf-1]
	fsym1 = ''
	for i in range(len(fseq1)):
		fsym1 += '#'
	begf = info[f'ch.{name}-wb']['exons'][-1][1]
	fseq2 = seq[begf:]
	fsym2 = ''
	for i in range(len(fseq2)):
		fsym2 += '#'
	esyms = []
	for e in info[f'ch.{name}-wb']['exons']:
		eseq = seq[e[0]-1:e[1]]
		esym = ''
		for i in range(len(eseq)): esym += '*'
		esyms.append(esym)
	isyms = []
	for i in info[f'ch.{name}-wb']['introns']:
		iseq = seq[i[0]-1:i[1]]
		isym = ''
		for i in range(len(iseq)): isym += '-'
		isyms.append(isym)
		
	sym_seq_wb = ''
	sym_seq_wb += fsym1
	for i in range(len(esyms)):
		sym_seq_wb += esyms[i]
		if i > len(isyms)-1: continue
		sym_seq_wb += isyms[i]
	sym_seq_wb += fsym2

	return sym_seq_wb

def make_apc_sym(iso, info, seq):

	first = info[iso]['exons'][0]
	if first[0] < first[1]:
		endf = first[0]
		fseq1 = seq[:endf-1]
		fsym1 = ''
		for i in range(len(fseq1)): fsym1 += '#'	
	if first[0] > first[1]:
		endf = max(first)
		fseq1 = seq[:endf-1]
		fsym1 = ''
		for i in range(len(fseq1)): fsym1 += '!'

	last = info[iso]['exons'][-1]
	if last[0] < last[1]:
		begf = last[1]
		fseq2 = seq[begf:]
		fsym2 = ''
		for i in range(len(fseq2)): fsym2 += '#'
	if last[0] > last[1]:
		begf = max(last)
		fseq2 = seq[begf-1:]
		fsym2 = ''
		for i in range(len(fseq2)): fsym2 += '!'
			
	esyms = []
	for e in info[iso]['exons']:
		if e[0] > e[1]:
			esyms.append('')
		if e[0] < e[1]:
			eseq = seq[e[0]-1:e[1]]
			esym = ''
			for i in range(len(eseq)): esym += '*'
			esyms.append(esym)

	isyms = []
	for i in info[iso]['introns']:
		iseq = seq[i[0]-1:i[1]]
		isym = ''
		for i in range(len(iseq)): isym += '-'
		isyms.append(isym)
		
	sym_seq_apc = ''
	sym_seq_apc += fsym1
	for i in range(len(esyms)):
		sym_seq_apc += esyms[i]
		if i > len(isyms)-1: continue
		sym_seq_apc += isyms[i]
	sym_seq_apc += fsym2

	return sym_seq_apc

def make_frame_sym(sym_seq_wb):

	c = 1
	frame = []
	for s in sym_seq_wb:
		if s == '#':
			frame.append('#')
		if s == '*':
			frame.append(c)
			c += 1
		if s == '-':
			frame.append('-')

	frame2 = ''
	for s in frame:
		if type(s) == int:
			if s%3 == 1: frame2 += '1'
			if s%3 == 2: frame2 += '2'
			if s%3 == 0: frame2 += '3'
		else:
			frame2 += s
	
	return frame2

def string_hyphs(seq, string):

	hyphs = ['-' for i in range(len(seq) - len(string) + 3)]
	line = string + ''.join(hyphs)

	return line
