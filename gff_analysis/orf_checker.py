import sys

fasta = sys.argv[1]
wb_gff = sys.argv[2]
apc_gff = sys.argv[3]

with open(fasta, 'r') as fp:
	seq = ''
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			gID = line
		else:
			seq += line

with open(wb_gff, 'r') as fp:
	count = 0
	for line in fp.readlines():
		line = line.rstrip()
		sline = line.split('\t')
		if sline[2] == 'CDS' and len(sline) == 9 and count == 0:
			print(sline)
			start = int(sline[3]) - 1
			print(seq[start:start+10])
			count += 1
		
print('#####')

with open(apc_gff, 'r') as fp:
	count = 0
	for line in fp.readlines():
		line = line.rstrip()
		sline = line.split('\t')
		if len(sline) == 9 and sline[2] == 'exon' and count == 0:
			print(sline[2])
			count += 1
print(seq[101:104])
print(seq[104-1:104+2])

# ATG #of bases STOP	
# stop codons: TAG, TAA, TGA

# organize different isoforms into dictionaries?
print('#####')

a_isos = {}
with open(apc_gff, 'r') as fp:
	count = 0
	for line in fp.readlines():
		line = line.rstrip()
		sline = line.split('\t')
		if len(sline) < 9: continue
		if sline[2] == 'gene': continue
		if sline[2] == 'mRNA':
			count += 1
			a_isos[count] = [sline]
		else:
			a_isos[count] += [sline]

for i in a_isos:
	print(i)
	for j in a_isos[i]:
		print(j)
		
