import sys
import modelib as ml
import isoform_fixed as isof

fp = sys.argv[1]

seq = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'

ppm, pwm = ml.read_pwm(fp)

dons, accs = ml.get_gtag(seq)

apc_isoforms, trials = ml.apc(dons, accs, 100, 3, 4, 5, seq)

for iso in apc_isoforms:
	break
print(iso)

d_seqs, a_seqs = ml.get_donacc_seqs(iso, seq)

print(pwm)
print(d_seqs)
print(a_seqs)

score_total = ml.get_pwm_score(d_seqs, pwm)

print(score_total)
print('*****')

for d_seq in d_seqs:
	print(d_seq)
	for i in range(len(d_seq)):
		print(i)
		print(pwm[i])
		'''
		print(nt)
		if nt == 'A': print(pwm[0][0])
		if nt == 'G': print(pwm[0][1])
		'''
