import sys
import modelib as ml


exons = ml.read_txt_seqs(sys.argv[1])

scores, probs = ml.make_mm(exons, 1)

print(probs)
print('*****')
print(scores)
