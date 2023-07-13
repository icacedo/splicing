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
		
ex_mm_pb, ex_mm_sc = read_exin_mm(exon_mm_tsv)

print(apc_isoforms)
print('*****')

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

#def get_mm_score(exin_seqs, exin_mm_model):
one_iso = apc_isoforms[0]
ex_seqs, in_seqs = get_exin_seqs(one_iso, seq)

one_seq = ex_seqs[0]

print(one_iso)
print(one_seq)
print('***')

k = 2
for i in range(len(one_seq)):
	if len(one_seq[i:i+k]) == k:
		kmer = one_seq[i:i+k]
		score = ex_mm_sc[kmer]
		print(kmer, score)

print(ex_mm_sc)

