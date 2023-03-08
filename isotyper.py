# use ch.49 as test gene
# intron 1: best predicted donor site out of frame (2k, 300k)
# intron 2 & 3: same don/acc site

import seqlib
import sys

fasta = sys.argv[1]
apc_info = 

for i,j in seqlib.read_fasta(fasta):
    print(i,j)


    

