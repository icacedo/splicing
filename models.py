import sys
import gzip
import modelib as ml

fp=gzip.open(sys.argv[1])
seqs=[]
for line in fp.readlines():
	# convert bytes to string
	line=line.decode('UTF=8')
	line=line.rstrip()
	seqs.append(line)

print(ml.make_pwm(seqs,6,len(seqs)))


	
