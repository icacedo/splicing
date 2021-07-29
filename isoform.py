import itertools
import sys

def short_intron(dons, accs, minintron):
	for d, a in zip(dons, accs):
		intron_length = a - d
		if intron_length < minintron: return True
	return False
	
def short_exon(dons, accs, minexon):
	for i in range(1, len(dons)):
		exon_beg = accs[i-1] + 1
		exon_end = dons[i] - 1
		exon_len = exon_end - exon_beg
		if exon_len < minexon: return True
	return False

"""
def all_probable_pwm(seq, minin, minex, dpwm, apwm, dt, at):
	dons = []
	accs = []
	for i in range(minex, len(seq) -minex):
		don = seq[i:i+len(dpwm)]
		acc = seq[i:i+len(apwm)]
		if pwmscore(don, dpwm) > dt: dons.append(i)
		if pwmscore(acc, apwm) > at: accs.append(i)	
	return generate_isoforms(dons, accs, minin, minex)
"""

def all_possible(seq, minin, minex):
	dons = []
	accs = []
	for i in range(minex, len(seq) -minex):
		if seq[i:i+2] == 'GT': dons.append(i)
		if seq[i:i+2] == 'AG': accs.append(i+1)
	return generate_isoforms(dons, accs, minin, minex)

def generate_isoforms(dons, accs, minin, minex):
	info = {
		'trials' : 0,
		'donors': len(dons),
		'acceptors': len(accs),
		'short_intron': 0,
		'short_exon': 0,
	}
	
	isoforms = []
	sites = min(len(dons), len(accs))
	for n in range(1, sites):
		for dsites in itertools.combinations(dons, n):
			for asites in itertools.combinations(accs, n):
				info['trials'] += 1
				
				# sanity checks
				if short_intron(dsites, asites, minin):
					info['short_intron'] += 1
					continue
				
				if short_exon(dsites, asites, minex):
					info['short_exon'] += 1
					continue
				
				# create isoform and save
				tx = []
				for d, a in zip(dsites, asites):
					tx.append({'beg':d, 'end':a})
				isoforms.append(tx)

	return isoforms, info
