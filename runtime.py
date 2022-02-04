import time
import seqlib
import sys
import os

for name, seq in seqlib.read_fasta(sys.argv[1]):
	with open('temp.fa', 'w') as fp:
		fp.write('>')
		fp.write(name)
		fp.write('\n')
		fp.write(seq)
	t0 = time.time()
	os.system('./geniso temp.fa | head -7 > temp.out')	
	t1 = time.time()
	print(t1-t0)
	
	with open('temp.out') as fp:
		lines = fp.readlines()
		name = lines[0]
