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


#np.savetxt(sys.argv[2],ml.make_pwm(seqs,6,len(seqs)))

# re-write make_simple_pwm to output in ACGT order
for i in ml.make_simple_pwm(seqs,6,len(seqs)):
	print(i)

'''
file=open(sys.argv[2],'w')
file.write(str(content))
file.close()
'''

	
