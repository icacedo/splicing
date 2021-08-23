import sys
import gzip

fp=gzip.open(sys.argv[1])
for line in fp.readlines():
	# convert bytes to string
	line=line.decode('UTF=8')
	line=line.rstrip()
	print(line)
	
