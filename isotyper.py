# use ch.49 as test gene
# intron 1: best predicted donor site out of frame (2k, 300k)
# intron 2 & 3: same don/acc site

import seqlib
import sys
import os
from grimoire.genome import Reader

# isoform structure info comes from the apc .gff3
# apc generates a gff file for the gene
# geniso makes a gff annotation for each isoform
fasta = sys.argv[1]
#apc_info = 

#for i,j in seqlib.read_fasta(fasta):
#    print(i,j)

