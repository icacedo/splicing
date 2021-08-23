import sys
import gzip
import modelib as ml
import numpy as np

fp=gzip.open(sys.argv[1])
seqs=[]
for line in fp.readlines():
	# convert bytes to string
	line=line.decode('UTF=8')
	line=line.rstrip()
	seqs.append(line)


np.savetxt(sys.argv[2],ml.make_pwm(seqs,6,len(seqs)))


	
