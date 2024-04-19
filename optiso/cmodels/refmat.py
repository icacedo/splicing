import sys

file = sys.argv[1]

count = 0
if 'len' in file:
	with open(file, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('%'): continue
			line = line.split('\t')
			print(line[0])
			count += 1
print(count)

if 'mm' in file:
	with open(file, 'r') as fp:
		count = 0
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('%'): continue
			line = line.split('\t')
			print(f'{line[0]} {line[1]}')
			count += 1
			if count%4 == 0: print(' ')

