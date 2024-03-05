# reformat length model arch/data/exon.len to .tsv
# reformat to read into score_all.py

import sys
import isoform_fixed as isof
import csv

fp = sys.argv[1]
filename = sys.argv[2]

len_probs = []
with open(fp, 'r') as lentxt:
	for line in lentxt.readlines():
		line = line.rstrip()
		if line.startswith('%'): continue
		len_probs.append(line)

len_model = isof.read_len(fp)
len_scores = len_model['val']

with open(filename, 'w', newline='') as tsvfile:
	writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
	writer.writerow(['% len P', 'len log2(P/expect)'])
	for i in range(len(len_probs)):
		writer.writerow([len_probs[i], len_scores[i]])
tsvfile.close()



















'''
#fp2 = sys.argv[2]
#fp3 = sys.argv[3]

# probabilities are converted to scores with these functions
thing = isof.read_len(fp)
print(thing['val'], thing['size'])
elen_scores = thing['val']

# prob to score using prob2score function
#thingy = isof.read_markov(fp2)
#print(thingy)

#thingay = isof.read_pwm(fp3)
#print(thingay)
'''








