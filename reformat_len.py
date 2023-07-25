# reformat length model arch/data/exon.len to .tsv
# reformat to read into score_all.py

import sys
import isoform_fixed as isof

fp = sys.argv[1]
'''
with open(fp, 'r') as lentxt:
	for line in lentxt.readlines():
		line = line.rstrip()
		print(line)
'''

fp2 = sys.argv[2]
fp3 = sys.argv[3]

thing = isof.read_len(fp)
print(thing['val'], thing['size'])

thingy = isof.read_markov(fp2)
print(thingy)

thingay = isof.read_pwm(fp3)
print(thingay)
