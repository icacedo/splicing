import isoform_fixed as isof
import sys
import copy
import modelib as ml

gff = sys.argv[1]

introns = isof.get_introns(gff)
for i in introns:
	print(i)

i1 = copy.deepcopy(introns)
print('#####')
for i in i1:
	print(i)


print('#####')

fasta = sys.argv[2]

seqid = None
seq = None
for seqid, seq in ml.read_fastas(fasta):
	seqid = seqid
	seq = seq

print(seq[100])
print(seq[726])





