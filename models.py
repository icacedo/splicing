# staring with the intron length model
# use arch/data/intron.txt.gz
# as an exercise
# try different smoothing techniques
# rectangular and with a slope
# or parabolic (goes to 0)
# start with linear vs rectangular
# rec, parameter is the length n of the rectangle
# linear, parameter is the slope m
# can be done in frequencies or sizes
# slope can create fractions
# can create a slope that goes by a set distance
# slope m, would be dependent on distance n
# how much distance? don't want to say intron of 0 is possible

import sys
import gzip

fp = sys.argv[1]

if fp.endswith(".gz"):
	with gzip.open(fp, 'r') as intfile:                       
	# everything is a byte string
	# can still slice and index like normal	
		for line in fp.readlines():
			line = line.rstrip()
			if isinstance(line, bytes):
				line = line.decode()
			print(line)

else:
	with open(fp, 'r') as intfile:
		for line in ifile.readlines():
			line = line.rstrip()
			print(line)
