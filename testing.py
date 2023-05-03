import modelib as ml
import sys
import isoform

'''
fp = sys.argv[1]

intbins = ml.get_intbins(fp)[0]

print(ml.rec_smoo(intbins))

print(ml.tri_smoo(intbins))
'''
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

'''
defaults:
maxs 100
minin 25
minex 25
flank 100
'''

minin = 3
minex = 4
flank = 5
maxs = 100

dons, accs = ml.get_gtag(seq1)
print(dons)
print(accs)
ml.apc(dons, accs, maxs, minin, minex, flank, seq1)

fp = sys.argv[1]
name, seq = next(isoform.read_fasta(fp))
print(name, seq)

