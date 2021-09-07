import sys
import gzip
import modelib as ml
import numpy as np

'''
fp = gzip.open(sys.argv[1])
seqs = []
for line in fp.readlines():
	# convert bytes to string
	line = line.decode('UTF=8')
	line = line.rstrip()
	seqs.append(line)
print(seqs)
'''

'''
# use this to create an output file for the numpy array
# can also use the direct output of ml.make_pwm() (no file created)
# this avoids having to recreate the numpy array from a file
np.savetxt(sys.argv[2],ml.make_pwm(seqs,6,len(seqs),0.001))
'''

'''
# read in pwm file, put back into array
fp = open("models_out/donor.pwm")
rows = []
arr = np.array
for line in fp.readlines():
	weight = ''
	for character in line:
		if character == ' ':
			rows.append(weight)
			weight = ''
		elif character == '\n':
			print(rows)
			rows = []
		else:
			weight += character
'''
arr = np.array([['4',5,6], [1,2,3]])
print(arr)







	
