import sys
import modelib as ml

# use ch.9940.fa for testing
fp = sys.argv[1]

seqid = None
seq = None
for seqid, seq in ml.read_fastas(fp):
	seqid = seqid
	seq = seq

#       0        9      16    22       31     38       47  
#                D      A     D        A      A          
seq1 = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
# iso1: AACATGACC        TTACCGTCACATTAGTTCGGAGCCCTATATA
# iso2: AACATGACC                       TTCGGAGCCCTATATA
# iso3: AACATGACC                              CCCTATATA
# iso4: AACATGACCGTTGCGAGTTACC          TTCGGAGCCCTATATA
# iso5: AACATGACCGTTGCGAGTTACC                 CCCTATATA
# iso6: AACATGACC        TTACC          TTCGGAGCCCTATATA
# iso7: AACATGACC        TTACC                 CCCTATATA

# defaults:
maxs = 100
minin = 25
minex = 25
flank = 100

# for seq1:
#maxs = 100
#minin = 3
#minex = 4
#flank = 5


dons, accs = ml.get_gtag(seq)
ml.apc(dons, accs, 100, 


