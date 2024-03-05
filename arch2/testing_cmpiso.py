import isoform_fixed as isof
import sys
import copy
import modelib as ml

gff1 = sys.argv[1]
gff2 = sys.argv[2]

introns1 = isof.get_introns(gff1)
introns2 = isof.get_introns(gff2)

print(introns1)
print(introns2)


for i in introns1:
	print(i)
print('#####')
for i in introns2:
	print(i)

print('****')

i1 = copy.deepcopy(introns1)
i2 = copy.deepcopy(introns2)

for i in i1:
	print(i)
print('*********')
for i in i2:
	print(i)


dist, details = isof.expdiff(introns1, introns2)

print(details)
print(dist)

for exon, p1, p2 in details:
	print(f'{exon[0]}\t{exon[1]}\t{p1:.6f}\t{p2:.6f}')

'''
i1 = copy.deepcopy(introns)
print('#####')
for i in i1:
	print(i)


print('#####')

fasta = sys.argv[3]

seqid = None
seq = None
for seqid, seq in ml.read_fastas(fasta):
	seqid = seqid
	seq = seq

print(seq[100])
print(seq[726])
'''




