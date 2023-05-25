import modelib as ml
import sys

seqs = ml.read_txt_seqs(sys.argv[1])

pfm = [{}for x in range(len(seqs[0]))]
print(pfm)

for s in seqs:
	for i in range(len(s)):
		print(s[i])

