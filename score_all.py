import sys
import modelib as ml

# ch.9940.fa
fp = sys.argv[1]

seqid = None
seq = None
for seqid, seq in ml.read_fastas(fp):
	seqid = seqid
	seq = seq
print(seqid)
print(seq)
