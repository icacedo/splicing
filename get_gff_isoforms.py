import sys

gff = sys.argv[1]

with open(gff) as fp:
	for line in fp.readlines():
		line = line.rstrip()
		print(line)




