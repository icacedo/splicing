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

def all_possible(seq, minin, minex):
	dons = []
	accs = []
	for i in range(minex, len(seq) -minex):
		if seq[i:i+2] == 'GT': dons.append(i)
		if seq[i:i+2] == 'AG': accs.append(i+1)

	info = {
		'trials' : 0,
		'donors': len(dons),
		'acceptors': len(accs),
		'no_intron' : 0,
		'unequal_count' : 0,
		'short_intron': 0,
		'short_exon': 0,
		'redundant': 0,
	}
	
	isoforms = []
	for n in range(0, len(dons)):
		for dsites in itertools.combinations(dons, n):
			for m in range(0, len(accs)):
				for asites in itertools.combinations(accs, m):
					info['trials'] += 1
					
					# sanity checks
					if len(dsites) == 0 or len(asites) == 0:
						info['no_intron'] += 1
						continue
					
					if len(dsites) != len(asites):
						info['unequal_count'] += 1
						continue
					
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

	# redundancy check -- shouldn't be any
	check = {}
	for iso in isoforms:
		shape = f'{iso}'
		if shape not in check: check[shape] = True
		else: info['redundant'] += 1

	return isoforms, info
