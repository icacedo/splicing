import sys

f = sys.argv[1]
fp = open(f)

count = 0
for line in fp.readlines():
	line = line.split('\t')
	print(line)
	count += 1
	if count == 10: sys.exit()
