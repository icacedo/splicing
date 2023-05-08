import sys

# get gt/ag sites from gff file
# use ch.9940

seq = ''
with open(sys.argv[2]) as faf:
	while True:
		line = faf.readline()
		line = line.rstrip()
		if not line: break
		print(line)
		if line.startswith('>'): continue
		seq += line
		
faf.close()		

print(seq)

with open(sys.argv[1]) as gff:
	while True:
		line = gff.readline()
		line = line.rstrip()
		if not line: break
		fields = line.split('\t')
		if fields[2] == 'intron':
			print(line)
			beg = int(fields[3])
			end = int(fields[4])
			print(beg, end)
			print(seq[beg-1:end])
			print(seq[beg-5:end+5])

gff.close()



