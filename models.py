import sys
import gzip
import itertools
import modelib as ml
import numpy as np
import isoform as iso
import seqlib as sql

for label, seq in sql.read_fasta(sys.argv[1]):
	print(label, seq)



#### re-write of the api algorithm, with more parameters ######################
# avg exon size in C. elegans: 200.7bp
# median: 123bp, smalest: 7bp, largest 7569bp
# avg intron size: 47bp
# median: 65bp, only 3 <= 25bp
# largest introns: 100912bp and 21230bp
# use GU/AG or GT/AG rule for donor/acceptor sites
# are the sequences in data/sequences/ translated mRNA or cDNA?
for i in itertools.combinations([1,2,3,4],2):
	print(i)
























	
