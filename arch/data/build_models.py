
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

exons = get_seqs('exon.txt.gz')
introns = get_seqs('intron.txt.gz')
accs = get_seqs('acceptor.txt.gz')
dons = get_seqs('donor.txt.gz')

# Exon Lengths
elen = isoform.create_len(exons, 15, 500)
isoform.write_len('exon.len', elen)

# Intron Lengths
ilen = isoform.create_len(introns, 5, 500)
isoform.write_len('intron.len', ilen)

# Acceptor PWM
apwm = isoform.create_pwm(accs)
isoform.write_pwm('acceptor.pwm', apwm)

# Donor PWM
dpwm = isoform.create_pwm(dons)
isoform.write_pwm('donor.pwm', dpwm)

# Exon Markov model
emm = isoform.create_markov(exons, 3, 0, 0)
isoform.write_markov('exon.mm', emm)

# Intron Markov model
imm = isoform.create_markov(introns, 3, 5, 6)
isoform.write_markov('intron.mm', imm)



