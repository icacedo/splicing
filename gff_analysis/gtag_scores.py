import modelib as ml
import sys

donor_pwm = sys.argv[1]
acceptor_pwm = sys.argv[2]

re_dppm, re_dpwm = ml.read_pwm(donor_pwm)
re_appm, re_apwm = ml.read_pwm(acceptor_pwm)

seq = 'AAAAGTGTCCAAAATATAAAGTACCACATCATTTCGAGCTCTCAGCGCGGC'
print(seq)
sites = []
dscores = []
for i in range(len(seq)):
	if seq[i:i+2] == 'GT':
		sites.append(i)
		dseq = seq[i:i+5]
		dpwm_score = ml.get_donacc_pwm_score(dseq, re_dpwm)
		dpwm_score = float('{:.5e}'.format(dpwm_score))
		dscores.append(dpwm_score)
	if seq[i:i+2] == 'AG':
		aseq = seq[i-4:i+2]
		apwm_score = ml.get_donacc_pwm_score(aseq, re_apwm)


for i in sites:
	print(i)
	for j in range(i+1):
		print(j)
	break


print(dscores[0])
for n in str(dscores[0]):
	print(n)







