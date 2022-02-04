import time
import seqlib
import sys
import os
'''
To filter out the sequences that too long to run using the APC algorithm 

To create a dataset with quickrunning data files
'''
fh = open('summary.tsv', 'w')

for name, seq in seqlib.read_fasta(sys.argv[1]):
	with open('temp.fa', 'w') as fp:
		fp.write('>')
		fp.write(name)
		fp.write('\n')
		fp.write(seq)
	t0 = time.time()
	os.system('./geniso temp.fa | head -7 > temp.out')	
	t1 = time.time()
	t = t1 - t0
	
	with open('temp.out') as fp:
		lines = fp.readlines()
		name = lines[0].split()[2]
		length = lines[1].split()[2]
		isoforms = lines[5].split()[2]
	
	fh.write(f'{name}	')
	fh.write(f'{length}	')
	fh.write(f'{isoforms}	')
	fh.write(f'{str(t)}	\n')
	
fh.close()
