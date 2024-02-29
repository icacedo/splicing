import modelib as ml
import sys

donor_pwm = sys.argv[1]
acceptor_pwm = sys.argv[2]

re_dppm, re_dpwm = ml.read_pwm(donor_pwm)
re_appm, re_apwm = ml.read_pwm(acceptor_pwm)

seq = 'AAAAGTGTCCAAAATATAAAGTACCACATCATTTCGAGCTCTCAGCGCGGC'

dscores = {}
ascores = {}
for i in range(len(seq)):
	if seq[i:i+2] == 'GT':
		dseq = seq[i:i+5]
		dpwm_score = ml.get_donacc_pwm_score(dseq, re_dpwm)
		dpwm_score = float('{:.5e}'.format(dpwm_score))
		dscores[i] = dpwm_score
	if seq[i:i+2] == 'AG':
		aseq = seq[i-4:i+2]
		apwm_score = ml.get_donacc_pwm_score(aseq, re_apwm)
		apwm_score = float('{:.5e}'.format(apwm_score))
		ascores[i] = apwm_score

print(dscores)
print(ascores)

lines = []
for site in dscores:
	line = ''
	for i in range(len(site)):
		line += '-'












