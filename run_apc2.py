import sys
import modelib as ml

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
# defaults for test seq:
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

dopwm_tsv_file = sys.argv[1]
acpwm_tsv_file = sys.argv[2]

donor_ppm, donor_pwm = read_pwm(dopwm_tsv_file)
acceptor_ppm, acceptor_pwm = read_pwm(acpwm_tsv_file)

for i in acceptor_ppm:
	print(i)
print('*****')
for i in acceptor_pwm:
	print(i)
print('*****')
for i in donor_ppm:
	print(i)
print('*****')
for i in donor_pwm:
	print(i)












