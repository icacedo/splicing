import sys
import gzip
import modelib as ml
import numpy as np
import isoform as iso
import allpossible_v2 as allv2

'''
fp = gzip.open(sys.argv[1])
seqs = []
for line in fp.readlines():
	# convert bytes to string
	line = line.decode('UTF=8')
	line = line.rstrip()
	seqs.append(line)
'''

'''
# use this to create an output file for the numpy array
np.savetxt(sys.argv[2],ml.make_pwm(seqs,6,len(seqs),0.001))
'''

'''
pwm_arr = ml.make_pwm(seqs,6,len(seqs),0.001)

site = 'GTACGC'
print(ml.site_score(pwm_arr, site))
'''

sites = ml.gff_reader(sys.argv[1], sys.argv[2], gtag=False)
don_sites = sites[0]
acc_sites = sites[1]
# Next step: all possible using gff sites

# compare output of allpossible_v2.py to isoform.py
seq = [seq for id, seq in iso.read_fasta(sys.argv[2])]
seq = seq[0]

print(iso.all_possible(seq, 10, 10, 10, 10, gff=sys.argv[1]))

print(allv2.all_possible(don_sites, acc_sites, 10, 10))


























	
