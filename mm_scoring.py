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

exon_mm_tsv = sys.argv[1]

def read_exin_mm(exin_mm_tsv):

	with open(exon_mm_tsv, 'r') as fp:
		re_mm_pb = {}
		re_mm_sc = {}
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			re_mm_pb[line[0]] = line[1]
			re_mm_sc[line[0]] = line[2]
		return re_mm_pb, re_mm_sc

def get_exin_seqs(isoform, seq):

	ex_seqs = []
	for exon in isoform['exons']:
		ex_beg = exon[0] 
		ex_end = exon[1] + 1
		exon_seq = seq[ex_beg:ex_end]
		ex_seqs.append(exon_seq)

	in_seqs = []
	for intron in isoform['introns']:
		in_beg = intron[0]
		in_end = intron[1] + 1
		intron_seq = seq[in_beg:in_end]
		in_seqs.append(intron_seq)

	return ex_seqs, in_seqs

def get_mm_score(exin_seqs, exin_mm):

	k = 0
	for key in exin_mm:
		string = key.split(' ')[2]
		k = int(string[0])
		break

	exin_score_total = 0
	for exin_seq in exin_seqs:
		exin_score = 0
		for i in range(len(exin_seq)):
			if len(exin_seq[i:i+k]) == k:
				kmer = exin_seq[i:i+k]
				score = exin_mm[kmer]
				exin_score += float(score)
		exin_score_total += exin_score
	return exin_score_total

ex_mm_pb, ex_mm_sc = read_exin_mm(exon_mm_tsv)
in_mm_sc = ex_mm_pb

for iso in apc_isoforms:
	ex_seqs, in_seqs = get_exin_seqs(iso, seq)
	emm_sc = get_mm_score(ex_seqs, ex_mm_sc)
	imm_sc = get_mm_score(in_seqs, in_mm_sc)
	iso['score'] = emm_sc + imm_sc
	print(iso)
	



