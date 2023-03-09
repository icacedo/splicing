# use ch.49 as test gene
# intron 1: best predicted donor site out of frame (2k, 300k)
# intron 2 & 3: same don/acc site

import seqlib
import sys
import os
from grimoire.genome import Reader

fasta = sys.argv[1]
#apc_info = 

#for i,j in seqlib.read_fasta(fasta):
#    print(i,j)

# where does intron information come from?
# copy pasted from apc2html.py

def get_isoforms(fasta,  params):
	os.system(f'isoformer {fasta} {params} > tmp.gff')
	genome = Reader(gff='tmp.gff', fasta=fasta)
	chrom = next(genome)
	for gene in chrom.ftable.build_genes():
		for tx in gene.transcripts():
			tx.use_longest_orf()
			yield tx

