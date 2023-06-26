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

# for seq1:
#maxs = 100
#minin = 3
#minex = 4
#flank = 5


dons, accs = ml.get_gtag(seq)
# ml.apc(dons, accs, maxs, minin, minex, flank, seq)
apc_isoforms, trials = ml.apc(dons, accs, 100, 3, 4, 5, seq)

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
		#print(seq[d_start:d_end], seq[a_start:a_end])
	print(d_seqs, a_seqs)
	break
	'''
	dstart = iso['introns'][0][0]
	dend = dstart + 5
	print(dstart, dend)
	print(seq[dstart:dend])
	break
	'''
print('***')
print(d_seqs, a_seqs)

donor_ppm = []
donor_pwm = []
count = 0
with open(sys.argv[2], 'r') as fp1:
	for line in fp1.readlines():
		line = line.rstrip()
		print(line.split('\t'))
		if count != 0 and  count <= 5:
			donor_ppm.append(line.split('\t'))
		if count >= 7:
			donor_pwm.append(line.split('\t'))
		count += 1
print('***')
print(donor_ppm)
print('***')
print(donor_pwm)
print('*****')

d_score = 0
for i in range(len(d_seqs)):
	print(d_seqs[i])
	for j in range(len(d_seqs[i])):
		if d_seqs[i][j] == 'A':
			print(d_seqs[i][j], i, j, donor_pwm[j][0])
		if d_seqs[i][j] == 'C':
			print(d_seqs[i][j], i, j, donor_pwm[j][1])
		if d_seqs[i][j] == 'G':
			print(d_seqs[i][j], i, j, donor_pwm[j][2])
		if d_seqs[i][j] == 'T':
			print(d_seqs[i][j], i, j, donor_pwm[j][3])
	break	











