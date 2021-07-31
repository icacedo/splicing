
import gzip
import isoform
import math

def get_seqs(filename):
	seqs = []
	with gzip.open(filename, 'rt') as fp:
		for seq in fp.readlines():
			seq = seq.rstrip()
			seqs.append(seq.rstrip())
	return seqs

exons = get_seqs('data/exon.txt.gz')

# Exon Lengths
ehist = isoform.create_hist(exons)
print(ehist)

# Intron Lengths
ihist = isoform.create_hist(introns)
print(ihist)

# Acceptor PWM
accs = get_seqs('data/acceptor.txt.gz')
apwm = isoform.create_pwm(accs)
print(apwm)

# Donor PWM
dons = get_seqs('data/donor.txt.gz')
dpwm = isoform.create_pwm(dons)
print(dpwm)

# Exon Markov model
emm = isoform.create_markov(exons, 3, 0, 0)
print(emm)

# Intron Markov model
introns = get_seqs('data/intron.txt.gz')
imm = isoform.create_markov(introns, 3, 5, 6)
print(imm)



