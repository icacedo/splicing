#!/usr/bin/env python3

import json
import os

from grimoire.genome import Reader
from grimoire.toolbox import translate_str

efp = open('exon.txt', 'w')
ifp = open('intron.txt', 'w')
dfp = open('donor.txt', 'w')
afp = open('acceptor.txt', 'w')

for region in os.listdir('favorites'): # from Lyman2020

	prefix = 'favorites/' + region + '/' + region
	f = open(prefix + '.json')
	meta = json.loads(f.read())
	f.close()

	genome = Reader(fasta=prefix+'.fa', gff=prefix+'.gff')
	chrom = genome.next()
	gene = chrom.ftable.build_genes()[0]
	
	tx = gene.transcripts()[0]
	
	for exon in tx.exons:
		efp.write(exon.seq_str())
		efp.write('\n')
	
	for intron in tx.introns:
		ifp.write(intron.seq_str())
		ifp.write('\n')
		dfp.write(intron.seq_str()[:5])   # GTAAG
		dfp.write('\n')
		afp.write(intron.seq_str()[-6:])  # TTTCAG
		afp.write('\n')
	
