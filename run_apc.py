import sys
import modelib as ml

# use ch.9940.fa for testing
fp = sys.argv[1]

seqid = None
seq = None
for seqid, seq in ml.read_fastas(fp):
	seqid = seqid
	seq = seq

#      0        9      16    22       31     38       47  
#               D      A     D        A      A          
seq = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'
# iso1: AACATGACC        TTACCGTCACATTAGTTCGGAGCCCTATATA
# iso2: AACATGACC                       TTCGGAGCCCTATATA
# iso3: AACATGACC                              CCCTATATA
# iso4: AACATGACCGTTGCGAGTTACC          TTCGGAGCCCTATATA
# iso5: AACATGACCGTTGCGAGTTACC                 CCCTATATA
# iso6: AACATGACC        TTACC          TTCGGAGCCCTATATA
# iso7: AACATGACC        TTACC                 CCCTATATA

# defaults:
maxs = 100
minin = 25
minex = 25
flank = 100

# for test seq:
maxs = 100
minin = 3
minex = 4
flank = 5

def read_pwm(ptsv):

	re_ppm = []
	re_pwm = []
	count = 0
	with open(ptsv, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			if count != 0 and  count <= 5:
				re_ppm.append(line.split('\t'))
			if count >= 7:
				re_pwm.append(line.split('\t'))
			count += 1
	return re_ppm, re_pwm

def get_pwm_score(da_seqs, da_pwm):

	for i in range(len(da_seqs)):
		da_score = 0
		for j in range(len(da_seqs[i])):
			if da_seqs[i][j] == 'A':
				da_score += float(da_pwm[j][0])
			if da_seqs[i][j] == 'C':
				da_score += float(da_pwm[j][1])
			if da_seqs[i][j] == 'G':
				da_score += float(da_pwm[j][2])
			if da_seqs[i][j] == 'T':
				da_score += float(da_pwm[j][3])
		yield da_seqs[i], da_score

dons, accs = ml.get_gtag(seq)
#apc_isos, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)
apc_isoforms, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)

donor_ppm, donor_pwm = read_pwm(sys.argv[2])

for iso in apc_isoforms:
	print(iso)
	d_seqs = []
	a_seqs = []
	for intron in iso['introns']:
		d_start = intron[0]
		d_end = d_start + 5
		a_end = intron[1] + 1
		a_start = a_end - 6
		d_seqs.append(seq[d_start:d_end])
		a_seqs.append(seq[a_start:a_end])
	print(d_seqs, a_seqs)
	for d_seq, d_score in get_pwm_score(d_seqs, donor_pwm):
		print(d_seq, d_score)


	












