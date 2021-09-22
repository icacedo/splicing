import sys
import gzip
import modelib as ml
import numpy as np
import isoform as iso

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
def gff_reader(gff):

	fp = open(gff)
	dons = []
	accs = []
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split()
		if line[2] == 'intron' and line[6] == '+': 
			if (int(line[3]) in dons) == False: 
				dons.append(int(line[3]))
			if (int(line[4]) in accs) == False: 
				accs.append(int(line[4]))
	
	return sorted(dons), sorted(accs)

print(gff_reader(sys.argv[1]))
print(iso.gff_sites(sys.argv[2], sys.argv[1]))


	
