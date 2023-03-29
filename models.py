# staring with the intron length model
# use arch/data/intron.txt.gz

import sys
import gzip

with gzip.open(sys.argv[1], 'r') as intfile:
	# everything is a byte string
	# can still slice and index like normal	
	for line in intfile.readlines():
		line = line.rstrip()
		if isinstance(line, bytes):
			line = line.decode()
		print(line)

