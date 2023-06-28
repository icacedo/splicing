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

exon_len_model = sys.argv[1]

count = 0
with open(exon_len_model, 'r') as fp:
	re_len_pdf = []
	re_len_sco = []
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('%'): continue
		if count <= 30:
			line = line.split('\t')
			print(line)
			re_len_pdf.append(line[0])
			re_len_sco.append(line[1])
			
		count += 1
print(re_len_pdf)
print(re_len_sco)

for iso in apc_isoforms:
	ex_lens = []
	print(iso)
	print(iso['exons'])
	for exon in iso['exons']:
		ex_len = exon[1] - exon[0] + 1
		print(exon, ex_len)
		ex_lens.append(ex_len)
	print(ex_lens)
	break
# double check if len model starts at 1 or 0
for elen in ex_lens:
	score = re_len_sco[elen]
	print(elen, score)



