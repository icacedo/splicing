import sys
import modelib as ml

# use ch.9940.fa for testing
'''
fp = sys.argv[1]

seqid = None
seq = None
for seqid, seq in ml.read_fastas(fp):
	seqid = seqid
	seq = seq
'''
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

def read_pwm(pwm_tsv):

	re_ppm = []
	re_pwm = []
	count = 0
	with open(pwm_tsv, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			if len(line.split('\t')) == 1:
				count +=1
				continue
			elif count == 1:
				re_ppm.append(line.split('\t'))
			elif count == 2:
				re_pwm.append(line.split('\t'))
	return re_ppm, re_pwm

def get_donacc_seqs(isoform):
	
	d_seqs = []
	a_seqs = []
	for intron in isoform['introns']:
		d_start = intron[0]
		d_end = d_start + 5
		a_end = intron[1] + 1
		a_start = a_end - 6
		d_seqs.append(seq[d_start:d_end])
		a_seqs.append(seq[a_start:a_end])
	return d_seqs, a_seqs

def get_pwm_score(da_seqs, da_pwm):

	da_score_total = 0
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
		da_score_total += da_score
		#print(da_seqs[i], da_score)
	return da_score_total

dons, accs = ml.get_gtag(seq)
apc_isoforms, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)

donor_pwm_tsv = sys.argv[1]
acceptor_pwm_tsv = sys.argv[2]

donor_ppm, donor_pwm = read_pwm(donor_pwm_tsv)
acceptor_ppm, acceptor_pwm = read_pwm(acceptor_pwm_tsv)

for iso in apc_isoforms:
	d_seqs, a_seqs = get_donacc_seqs(iso)
	dpwm_sc = get_pwm_score(d_seqs, donor_pwm)
	apwm_sc = get_pwm_score(a_seqs, acceptor_pwm)
	iso['score'] = dpwm_sc + apwm_sc
	print(iso)










