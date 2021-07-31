
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
			p = count[i][nt] / len(seqs)
			if p == 0:
				pwm[i][nt] = -99
			else:
				pwm[i][nt] = math.log2(p/0.25)
	return pwm

def read_pwm(file):
	# open file
	# read pwm
	# return pwm
	pass

def write_pwm(file, pwm):
	# open file for writing
	# write pwm
	pass

def score_pwm(seq, pwm):
	# assert seq is same length as pwm
	# return score
	pass

####################
## LENGTH SECTION ##
####################

def create_hist(seqs):
	hist = []
	h2 = {}
	for seq in seqs:
		n = len(seq)
		while len(hist) < n+1 :
			hist.append(0)
		if n not in h2: h2[n] = 0
		h2[n] += 1
		hist[n] += 1
	for i in range(len(hist)):
		v = None
		if i not in h2: v = '*'
		else: v = h2[i]
		print(i, hist[i], v)

def read_hist(file):
	# open file
	# read hist
	# return hist
	pass

def write_hist(file, hist):
	# open file for writing
	# write hist
	pass

def score_hist(seq, pwm):
	# return score
	pass

##########################
## MARKOV MODEL SECTION ##
##########################

def create_markov(seqs, order, beg, end):
	mm = {}
	for seq in seqs:
		for i in range(beg+order, len(seq) - end):
			ctx = seq[i-order:i]
			nt = seq[i]
			if ctx not in mm: mm[ctx] = {'A':0, 'C':0, 'G':0, 'T':0}
			mm[ctx][nt] += 1
	return mm

def read_markov(seqs):
	# open file
	# read model
	# return model
	pass

def write_markov(file, mm):
	# open file for writing
	# write model
	pass

def score_makov(seq, mm):
	# build score
	# return score
	pass

################################
## ISOFORM GENERATION SECTION ##
################################

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

def all_probable(seq, mini, mine, maxs, ignore,
		ihist=None, ehist=None, dpwm=None, apwm=None, imm=None, emm=None):
	# looks like all_possible but with optional filters
		# for acceptor and donor matches to pwms
		# for probabilistic lengths of introns and exons
		# for markov models of intron and exon composition
		# for final build?
	pass

def all_possible(seq, minin, minex, maxs, ignore):
	dons = []
	accs = []
	for i in range(ignore, len(seq) -ignore):
		if seq[i:i+2]   == 'GT': dons.append(i)
		if seq[i-1:i+1] == 'AG': accs.append(i)

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
				
				if short_exon(dsites, asites, minex):
					info['short_exon'] += 1
					continue
				
				# create isoform and save
				tx = []
				for d, a in zip(dsites, asites):
					tx.append({'beg':d, 'end':a})
				isoforms.append(tx)

	return isoforms, info

