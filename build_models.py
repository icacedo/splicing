
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
introns = get_seqs('data/intron.txt.gz')
accs = get_seqs('data/acceptor.txt.gz')
dons = get_seqs('data/donor.txt.gz')

# Exon Lengths
elen = isoform.create_len(exons, 25, 1000)
isoform.write_len('data/exon.len', elen)

# Intron Lengths
ilen = isoform.create_len(introns, 5, 1000)
isoform.write_len('data/intron.len', ilen)

# Acceptor PWM
apwm = isoform.create_pwm(accs)
isoform.write_pwm('data/acceptor.pwm', apwm)

# Donor PWM
dpwm = isoform.create_pwm(dons)
isoform.write_pwm('data/donor.pwm', dpwm)

# Exon Markov model
emm = isoform.create_markov(exons, 3, 0, 0)
isoform.write_markov('data/exon.mm', emm)

# Intron Markov model
imm = isoform.create_markov(introns, 3, 5, 6)
isoform.write_markov('data/intron.mm', imm)



