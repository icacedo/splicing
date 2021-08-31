
import itertools
import math
import random
import sys

#####################
## UTILITY SECTION ##
#####################

def randseq(n):
	seq = ''
	for i in range(n):
		seq += random.choice('ACGT')
	return seq

def get_filepointer(filename):
	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)
	return fp

def read_fasta(filename):

	name = None
	seqs = []

	fp = get_filepointer(filename)

	for line in fp.readlines(): # consider changing this
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()


#################
## PWM SECTION ##
#################

def create_pwm(seqs):
	count = []
	for seq in seqs:
		for i, nt in enumerate(seq):
			if len(count) <= i:
				count.append({'A':0, 'C': 0, 'G': 0, 'T': 0})
			count[i][nt] += 1

	pwm = [{} for i in range(len(count))]
	for i in range(len(count)):
		for nt in count[i]:
			pwm[i][nt] = count[i][nt] / len(seqs)
	return pwm

def write_pwm(file, pwm):
	with open(file, 'w') as fp:
		fp.write(f'# PWM {file} {len(pwm)}\n')
		for pos in pwm:
			for nt in pos:
				fp.write(f'{pos[nt]:.6f} ')
			fp.write('\n')

def read_pwm(file):
	# open file
	# read pwm
	# return pwm
	pass

def score_pwm(seq, pwm):
	# assert seq is same length as pwm
	# return score
	pass

####################
## LENGTH SECTION ##
####################

def create_len(seqs, floor, limit):
	count = []
	sum = 0
	icount = 0
	for seq in seqs:
		n = len(seq)
		sum += n
		icount += 1
		while len(count) < n+1 :
			count.append(0)

		count[n] += 1

	# rectangular smoothing
	r = 5 # 5 on each side
	smooth = [0 for i in range(len(count))]
	for i in range(r, len(count) -r):
		for j in range(-r, r+1):
			smooth[i+j] += count[i]

	for i in range(floor):
		smooth[i] = 0
	smooth = smooth[:limit]

	# model
	model = []
	total = 0
	for v in smooth: total += v
	for v in smooth: model.append(v/total)

	return model

def write_len(file, hist):
	with open(file, 'w') as fp:
		fp.write(f'# LEN {file} {len(hist)}\n')
		for val in hist:
			fp.write(f'{val:.6f}\n')

def read_len(file):
	# open file
	# read hist
	# return hist
	pass

def score_len(seq, pwm):
	# return score
	pass

##########################
## MARKOV MODEL SECTION ##
##########################

def create_markov(seqs, order, beg, end):
	count = {}
	for seq in seqs:
		for i in range(beg+order, len(seq) - end):
			ctx = seq[i-order:i]
			nt = seq[i]
			if ctx not in count: count[ctx] = {'A':0, 'C':0, 'G':0, 'T':0}
			count[ctx][nt] += 1

	# these need to be probabilities
	mm = {}
	for kmer in count:
		mm[kmer] = {}
		total = 0
		for nt in count[kmer]: total += count[kmer][nt]
		for nt in count[kmer]: mm[kmer][nt] = count[kmer][nt] / total

	return mm

def write_markov(file, mm):
	with open(file, 'w') as fp:
		fp.write(f'# MM {file} {len(mm)*4}\n')
		for kmer in sorted(mm):
			#fp.write(f'{kmer} ')
			for v in mm[kmer]:
				fp.write(f'{kmer}{v} {mm[kmer][v]:.6f}\n')
			fp.write('\n')

def read_markov(seqs):
	# open file
	# read model
	# return model
	pass

def score_makov(seq, mm):
	# build score
	# return score
	pass

################################
## TRANSCRIPT SCORING SECTION ##
################################

def score_apwm(tx, pwm):
	pass

def score_dpwm(tx, pwm):
	pass

def score_elen(tx, model):
	pass

def score_ilen(tx, model):
	pass

def score_emm(tx, mm):
	pass

def score_imm(tx, mm):
	pass

################################
## ISOFORM GENERATION SECTION ##
################################

def short_intron(dons, accs, min):
	for d, a in zip(dons, accs):
		intron_length = a - d + 1
		if intron_length < min: return True
	return False

def short_exon(dons, accs, seqlen, flank, min):

	# first exon
	exon_beg = flank + 1
	exon_end = dons[0] -1
	exon_len = exon_end - exon_beg + 1
	if exon_len < min: return True

	# last exon
	exon_beg = accs[-1] + 1
	exon_end = seqlen - flank + 1
	exon_len = exon_end - exon_beg + 1
	if exon_len < min: return True

	# interior exons
	for i in range(1, len(dons)):
		exon_beg = accs[i-1] + 1
		exon_end = dons[i] - 1
		exon_len = exon_end - exon_beg
		if exon_len < min: return True
	return False

def gtag_sites(seq, flank, minex):
	dons = []
	accs = []
	for i in range(flank + minex, len(seq) -flank -minex):
		if seq[i:i+2]   == 'GT': dons.append(i)
		if seq[i-1:i+1] == 'AG': accs.append(i)
	return dons, accs

def gff_sites(seq, gff):
	dond = {}
	accd = {}
	with open(gff) as fp:
		for line in fp.readlines():
			f = line.split()
			if f[2] != 'intron': continue
			dond[int(f[3]) -1] = True
			accd[int(f[4]) -1] = True

	dons = list(dond.keys())
	accs = list(accd.keys())
	dons.sort()
	accs.sort()

	return dons, accs

def build_mRNA(seq, beg, end, dons, accs):
	assert(beg <= end)
	assert(len(dons) == len(accs))

	tx = {
		'seq': seq,
		'beg': beg,
		'end': end,
		'exons': [],
		'introns': [],
		'score': 0
	}

	if len(dons) == 0:
		tx['exons'].append((beg, end))
		return tx

	# introns
	for a, b in zip(dons, accs):
		tx['introns'].append((a, b))

	# exons
	tx['exons'].append((beg, dons[0] -1))
	for i in range(1, len(dons)):
		a = accs[i-1] + 1
		b = dons[i] -1
		tx['exons'].append((a, b))
	tx['exons'].append((accs[-1] +1, end))

	return tx

def all_possible(seq, minin, minex, maxs, flank, dpwm=None, apwm=None,
		emm=None, imm=None, elen=None, ilen=None, gff=None):

	if gff: dons, accs = gff_sites(seq, gff)
	else:   dons, accs = gtag_sites(seq, flank, minex)

	info = {
		'trials' : 0,
		'donors': len(dons),
		'acceptors': len(accs),
		'short_intron': 0,
		'short_exon': 0,
	}

	isoforms = []
	sites = min(len(dons), len(accs), maxs)
	for n in range(1, sites+1):
		for dsites in itertools.combinations(dons, n):
			for asites in itertools.combinations(accs, n):
				info['trials'] += 1

				# sanity checks
				if short_intron(dsites, asites, minin):
					info['short_intron'] += 1
					continue

				if short_exon(dsites, asites, len(seq), flank, minex):
					info['short_exon'] += 1
					continue

				# create isoform and save
				score = 0
				tx = build_mRNA(seq, flank, len(seq) -flank, dsites, asites)
				if apwm: score += score_apwm(apwm, tx)
				if dpwm: score += score_dpwm(dpwm, tx)
				if elen: score += score_elen(elen, tx)
				if ilen: score += score_ilen(ilen, tx)
				if emm:  score += score_emm(emm, tx)
				if imm:  score += score_imm(imm, tx)
				tx['score'] = score
				isoforms.append(tx)

	return isoforms, info
