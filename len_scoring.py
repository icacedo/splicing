import sys
import modelib as ml

#      0        9      16    22       31     38       47  
#               D      A     D        A      A          
seq = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
# iso1: AACATGACC        TTACCGTCACATTAGTTCGGAGCCCTATATA
# iso2: AACATGACC                       TTCGGAGCCCTATATA
# iso3: AACATGACC                              CCCTATATA
# iso4: AACATGACCGTTGCGAGTTACC          TTCGGAGCCCTATATA
# iso5: AACATGACCGTTGCGAGTTACC                 CCCTATATA
# iso6: AACATGACC        TTACC          TTCGGAGCCCTATATA
# iso7: AACATGACC        TTACC                 CCCTATATA
# for test seq:
maxs = 100
minin = 3
minex = 4
flank = 5

dons, accs = ml.get_gtag(seq)
apc_isoforms, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)

exon_len_tsv = sys.argv[1]

def read_exin_len(exin_len_tsv):

	with open(exin_len_tsv, 'r') as fp:
		re_len_pdf = []
		re_len_sco = []
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('%'): continue
			line = line.split('\t')
			re_len_pdf.append(line[0])
			re_len_sco.append(line[1])
	return re_len_pdf, re_len_sco

def get_exin_lengths(isoform):
	
	ex_lens = []
	for exon in isoform['exons']:
		ex_len = exon[1] - exon[0] + 1
		ex_lens.append(ex_len)
	in_lens = []
	for intron in isoform['introns']:
		in_len = intron[1] - intron[0] +1
		in_lens.append(in_len)
	return ex_lens, in_lens

def get_len_score(exin_lens, exin_len_model):

	exin_score_total = 0
	for length in ex_lens:
		exin_score = exin_len_model[length]
		exin_score_total += float(exin_score)
	return(exin_score_total)

exon_len_pdf, exon_len_sco = read_exin_len(exon_len_tsv)

for iso in apc_isoforms:	
	ex_lens, in_lens = get_exin_lengths(iso)
	elen_sc = get_len_score(ex_lens, exon_len_sco)
	iso['score'] = elen_sc 
	print(iso)
	





