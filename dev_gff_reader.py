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

dons = []
accs = []
with open(sys.argv[1]) as gff:
	while True:
		line = gff.readline()
		line = line.rstrip()
		if not line: break
		fields = line.split('\t')
		if fields[2] == 'intron':
			print(line)
			beg = int(fields[3]) - 1
			end = int(fields[4]) - 1
			print(beg, end)
			dons.append(beg)
			accs.append(end)
			print(seq[beg:end+1])
			print(seq[beg:beg+2]) # get GT
			print(seq[end-1:end+1]) # get AG
	
gff.close()	

import modelib

print(dons, accs)
print('***')
dons1, accs1 = modelib.get_gtag(seq)
print(dons1, accs1)

# next test ian's code




