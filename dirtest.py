import os
import sys

apc_dir = sys.argv[1]

#print(os.listdir(apc_dir))
#print(apc_dir)

for file in os.listdir(apc_dir):
	print(file)
	fname = apc_dir + file
	print(fname)
	with open(fname, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			print(line)
	break
