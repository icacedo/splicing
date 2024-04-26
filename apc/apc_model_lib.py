import sys
import gzip
import math
import openturns as ot
from itertools import combinations

################################
##### File Reading Section #####
################################

def read_fastas(fastafile):

	with open(fastafile) as fp:

		seqid = ''
		seq = ''	

		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('>'): seqid = line[1:]
			else: seq += line
		yield (seqid, seq)			
	fp.close()		

def read_gff_sites(seq, gff, gtag=True):

	dons = []
	accs = []

	with open(gff) as fp:
		while True:
			line = fp.readline()
			line = line.rstrip()
			if not line: break
			fields = line.split('\t')
			if fields[2] == 'intron':
				beg = int(fields[3]) - 1
				end = int(fields[4]) - 1
				if gtag:
					if seq[beg:beg+2] == 'GT': 
						dons.append(beg)
					if seq[end-1:end+1] == 'AG': 
						accs.append(end)
	fp.close()	
	return sorted(set(dons)), sorted(set(accs))

# read exon/intron or donor/acceptor sequences file
def read_txt_seqs(fp):

	tn_seqs = [] # training sequences

	if fp.endswith(".gz"):
		with gzip.open(fp, 'r') as tnseqfile:
			tn_seq = ''
			for line in tnseqfile.readlines():
				line = line.rstrip()
				if isinstance(line, bytes):
					line = line.decode()
				tn_seq = line
				tn_seqs.append(tn_seq)
		tnseqfile.close()
		return tn_seqs 

	else:
		with open(fp, 'r') as tnseqfile:
			tn_seq = ''
			for line in tnseqfile.readlines():
				line = line.rstrip()
				tn_seq = line
				tn_seqs.append(tn_seq)
			tnseqfile.close()
		return tn_seqs

################################
##### Length Model Section #####
################################

# returns a count/frequency  histogram or raw counts
# uses a list of exons/introns as input
def get_exinbins(exinseqs, nbins=None, pre=None):

	lines = []
	total_obs = 0
	for l in exinseqs:
		lines.append(l)
		total_obs += 1

	sizes = []
	for exin in lines:
		sizes.append(len(exin))
	
	freqs = []
	total = len(sizes)
	for count in sizes:
		fq = count/total
		if pre:
			fq2 = f"{fq:.{pre}f}"
			freqs.append(float(fq2))
		else: freqs.append(float(fq))
	
	count_bins = [0 for x in range(max(sizes)+1)]
	for i in range(len(sizes)):
		count_bins[sizes[i]] += 1

	freq_bins = []
	for i in range(len(count_bins)):
		fqb = count_bins[i]/total_obs
		if pre:	
			fqb2 = f"{fqb:.{pre}f}"
			freq_bins.append(float(fqb2))	
		else: freq_bins.append(float(fqb))	

	# only append up to nbins, for testing
	return count_bins[:nbins], freq_bins[:nbins], sizes, freqs

##### smoothing #####

# definitions of rectangular and triangular smoothing found here:
# https://terpconnect.umd.edu/~toh/spectrum/Smoothing.html

# rectangular smoothing
def rec_smoo(intbins, m=5, pre=None):
	
	smoodata = []
	for i in range(len(intbins)):
		m2 = int((m/2) + 0.5 - 1)
		bef = intbins[i-m2:i]
		now = intbins[i]
		aft = intbins[i+1:i+m2+1]
		total = sum(bef) + now + sum(aft)
		smoopt = total/m
		if pre:
			fmat = f'{smoopt:.{pre}f}'
			smoodata.append(fmat)
		else:
			smoodata.append(smoopt)
		#print(bef, now, aft, total, smoopt)
	
	return smoodata

# triangular smoothing
def tri_smoo(intbins, m=5, pre=None):	

	smoodata = []
	m2 = int((m/2) + 0.5 - 1)
	for i in range(len(intbins)):
		inx = i-m2
		if inx < 0: inx = 0
		bef = intbins[inx:i]
		now = intbins[i]
		aft = intbins[i+1:i+m2+1]
		nowx = now * (m2 + 1)
		tbefx = 0
		cbefx = 0
		for j in range(len(bef)):
			coefb = m2 - 1 + j
			befx = bef[j] * coefb
			tbefx += befx
			cbefx += coefb
		taftx = 0
		caftx = 0
		for k in range(len(aft)):
			coefa = m2 - k
			aftx = aft[k] * coefa
			taftx += aftx
			caftx += coefa
		total_coef = (m2 + 1) + cbefx + caftx
		total = nowx + tbefx + taftx
		smoopt = total/total_coef
		if pre: 
			fmat = f'{smoopt:.{pre}f}'
			smoodata.append(fmat)
		else:
			smoodata.append(smoopt)
		#print(bef, now, aft, total, smoopt)

	return smoodata

##### curve fitting #####

def frechet_pdf(x, a, b, g):
	if x < g: return 0
	z = (x-g)/b
	term1 = (a/b)
	term2 = z**(-1-a)
	term3 = math.exp(-z**-a)
	return term1 * term2 * term3

def fdist_params(exinseqs, nbins=None, pre=None, size_limit=None):

	data = get_exinbins(exinseqs, nbins=None, pre=None)[2]
	
	if size_limit:
		sample = ot.Sample([[x] for x in data if x < size_limit])
	else:
		size_limit = max(data)
		sample = ot.Sample([[x] for x in data if x < size_limit])	

	distFrechet = ot.FrechetFactory().buildAsFrechet(sample)

	a = distFrechet.getAlpha()
	b = distFrechet.getBeta()
	g = distFrechet.getGamma()
	
	return data, a, b, g, size_limit

# at size_limit 500-max exon/intron fit frechet dist
# lower than 500, no longer frechet dist
# max intron size is 5862
# max exon size is 2921
def memoize_fdist(data, a, b, g, size_limit, pre=None):
		
	x_values = []
	y_values = []
	y_scores = []
	# is this necessary for extreme value distribution?
	expect = 1/size_limit
	for i in range(min(len(data), size_limit)):
		x_values.append(i)
		y = frechet_pdf(i, a, b, g)
		if pre:
			y2 = f"{y:.{pre}f}"
			y_values.append(y2)
			if y == 0: y_scores.append(-100)
			else:
				#ys = math.log2(y) 
				ys = math.log2(y/expect)
				ys2 =  f"{ys:.{pre}f}"
				y_scores.append(ys2)
		else: 
			y_values.append(y)
			if y == 0: y_scores.append(-100)
			#else: y_scores.append(math.log2(y/expect))
			else: y_scores.append(math.log2(y/expect))			
	
	data = {
		'x': x_values,
		'y': y_values
	}
	
	# only scores are useful?
	# y beore log transforming
	return y_scores, y_values

##### len model scoring #####

def read_exin_len(exin_len_tsv):

	with open(exin_len_tsv, 'r') as fp:
		re_len_pdf = []
		re_len_sco = []
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('%'): continue
			line = line.split('\t')
			re_len_pdf.append(line[0])
			re_len_sco.append(line[1])
	return re_len_pdf, re_len_sco

def read_len_params(exin_len_tsv):
	
	with open(exin_len_tsv, 'r') as fp:
		a = None
		b = None
		g = None
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('% EVD params:'):
				line = line.split(' ')
				a = float(line[4])
				b = float(line[6])
				g = float(line[8])
			else: break
		return a, b, g

def get_exin_len_score(exin, exin_len_model, a, b, g):

	exin_len = exin[1] - exin[0] + 1

	if exin_len < len(exin_len_model):
		exin_len_score = exin_len_model[exin_len]
	else:
		exin_prob = frechet_pdf(exin_len, a, b, g)
		expect = 1/exin_len
		exin_len_score = math.log2(exin_prob/expect)
	return float(exin_len_score)

# old method below
'''
def get_exin_lengths(isoform):
	
	ex_lens = []
	for exon in isoform['exons']:
		ex_len = exon[1] - exon[0] + 1
		ex_lens.append(ex_len)
	in_lens = []
	for intron in isoform['introns']:
		in_len = intron[1] - intron[0] + 1
		in_lens.append(in_len)
	return ex_lens, in_lens

def get_len_score(exin_lens, exin_len_model, a, b, g):

	exin_score_total = 0
	for length in exin_lens:
		if length < len(exin_len_model):
			exin_score = exin_len_model[length]
			exin_score_total += float(exin_score)
		else: 
			exin_prob = frechet_pdf(length, a, b, g)
			# not sure how to do the expectation for longer exons/introns
			expect = 1/length
			exin_score = math.log2(exin_prob/expect)
			exin_score_total += float(exin_score)
			
	return exin_score_total
'''
################################
##### Markov Model Section #####
################################

def make_mm(exinseqs, order=3):

	context = {}
	for seq in exinseqs:
		for i in range(len(seq)-order):
			prev = seq[i:i+order]
			now = seq[i+order]
			if prev not in context:
				context[prev] = now
			else:
				context[seq[i:i+order]] += now
	mm_scores = {}
	mm_probs = {}
	for nts in sorted(context):	
		A = 0
		C = 0
		G = 0
		T = 0
		for nt in context[nts]:
			if nt == 'A': A += 1
			if nt == 'C': C += 1
			if nt == 'G': G += 1
			if nt == 'T': T += 1
			d = int(len(context[nts]))
		if nts not in mm_scores:
			if A/d == 0: ad = -100
			else: ad = math.log2((A/d)/0.25)
			if C/d == 0: cd = -100
			else: cd = math.log2((C/d)/0.25)
			if G/d == 0: gd = -100 
			else: gd = math.log2((G/d)/0.25)
			if T/d == 0: td = -100
			else: td = math.log2((T/d)/0.25)	
			mm_scores[nts] = (ad, cd, gd, td)
		if nts not in mm_probs:
			mm_probs[nts] = (A/d, C/d, G/d, T/d)
	return mm_scores, mm_probs, order

##### Markov model scoring #####

def read_exin_mm(exin_mm_tsv):

	with open(exin_mm_tsv, 'r') as fp:
		re_mm_pb = {}
		re_mm_sc = {}
		for line in fp.readlines():
			if line.startswith('%'): continue
			line = line.rstrip()
			line = line.split('\t')
			re_mm_pb[line[0]] = line[1]
			re_mm_sc[line[0]] = line[2]
		return re_mm_pb, re_mm_sc

def get_exin_mm_score(exin, seq, exin_mm, dpwm=None, apwm=None):

	beg = exin[0]
	end = exin[1] + 1
	exin_seq = seq[beg:end]

	k = 0
	for key in exin_mm:
		k = len(key)
		break

	if dpwm and apwm:
		exin_seq = exin_seq[len(dpwm):-len(apwm)]
	
	exin_mm_score = 0
	for i in range(len(exin_seq)):
		if len(exin_seq[i:i+k]) == k:
			kmer = exin_seq[i:i+k]
			exin_mm_score += float(exin_mm[kmer])

	return float(exin_mm_score)
	
# old method below
'''
def get_exin_seqs(isoform, seq):

	ex_seqs = []
	for exon in isoform['exons']:
		ex_beg = exon[0] 
		ex_end = exon[1] + 1
		exon_seq = seq[ex_beg:ex_end]
		ex_seqs.append(exon_seq)

	in_seqs = []
	for intron in isoform['introns']:
		in_beg = intron[0]
		in_end = intron[1] + 1
		intron_seq = seq[in_beg:in_end]
		in_seqs.append(intron_seq)

	return ex_seqs, in_seqs

def get_mm_score(exin_seqs, exin_mm, dpwm=None, apwm=None):
	
	k = 0
	for key in exin_mm:
		k = len(key)
		break
	
	exin_seqs2 = []	
	if dpwm and apwm:
		for in_seq in exin_seqs:
			in_seq = in_seq[len(dpwm):-len(apwm)]
			exin_seqs2.append(in_seq)
	else:
		exin_seqs2 = exin_seqs

	exin_score_total = 0
	for exin_seq in exin_seqs2:
		exin_score = 0
		for i in range(len(exin_seq)):
			if len(exin_seq[i:i+k]) == k:
				kmer = exin_seq[i:i+k]
				score = exin_mm[kmer]
				exin_score += float(score)
		exin_score_total += exin_score
	return exin_score_total
'''
#######################
##### PWM section #####
#######################

def make_pwm(seqs):

	pfm = [{'A': 0, 'C': 0, 'G': 0, 'T':0} for x in range(len(seqs[0]))]
	for i in range(len(seqs)):
		for j in range(len(seqs[i])):
			pfm[j][seqs[i][j]] += 1

	ppm = [{'A': 0, 'C': 0, 'G': 0, 'T':0} for x in range(len(pfm))]
	for i in range(len(pfm)):
		for n in pfm[i]:
			frequency = pfm[i][n]/len(seqs)		
			ppm[i][n] = frequency	

	pwm = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for x in range(len(ppm))]
	for i in range(len(ppm)):	
		for n in ppm[i]:
			if ppm[i][n] == 0:
				pwm[i][n] = -100
			else:
				weight = math.log2(ppm[i][n]/0.25)		
				pwm[i][n] = weight
	
	return pwm, ppm

##### PWM scoring ###

def read_pwm(pwm_tsv):

	re_pwm = []
	re_ppm = []
	count = 0
	with open(pwm_tsv, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('%'): 
				count += 1
				continue	
			elif count == 1:
				re_pwm.append(line.split('\t'))
			elif count == 2:
				re_ppm.append(line.split('\t'))
	return re_ppm, re_pwm

def get_donacc_seq(intron, seq):

	d_start = intron[0]
	d_end = d_start + 5
	a_end = intron[1] + 1
	a_start = a_end -6
	d_seq = seq[d_start:d_end]
	a_seq = seq[a_start:a_end]
	
	return d_seq, a_seq

def get_donacc_pwm_score(donacc, pwm):
	
	da_score = 0
	count = 0
	for i in range(len(donacc)):
		if donacc[i] == 'A':
			da_score += float(pwm[count][0])
		if donacc[i] == 'C':
			da_score += float(pwm[count][1])
		if donacc[i] == 'G':
			da_score += float(pwm[count][2])
		if donacc[i] == 'T':
			da_score += float(pwm[count][3])
		count += 1

	return da_score

# old method
'''
def get_donacc_seqs(isoform, seq):
	
	d_seqs = []
	a_seqs = []
	for intron in isoform['introns']:
		d_start = intron[0]
		d_end = d_start + 5
		a_end = intron[1] + 1
		a_start = a_end - 6
		d_seqs.append(seq[d_start:d_end])
		a_seqs.append(seq[a_start:a_end])
	return d_seqs, a_seqs

def get_pwm_score(da_seqs, da_pwm):

	da_score_total = 0
	for i in range(len(da_seqs)):
		da_score = 0
		for j in range(len(da_seqs[i])):
			if da_seqs[i][j] == 'A':
				da_score += float(da_pwm[j][0])
			if da_seqs[i][j] == 'C':
				da_score += float(da_pwm[j][1])
			if da_seqs[i][j] == 'G':
				da_score += float(da_pwm[j][2])
			if da_seqs[i][j] == 'T':
				da_score += float(da_pwm[j][3])
		da_score_total += da_score
		#print(da_seqs[i], da_score)
	return da_score_total
'''
##### other #####

def get_entropy(probs):

	h = 0
	for p in probs:
		h -= p * math.log2(p)
	return h

#######################
##### APC Section #####
#######################

def get_gtag(seq):

	dons = []
	accs = []
	for i in range(len(seq)):
		if seq[i:i+2] == 'GT':
			dons.append(i)
		if seq[i:i+2] == 'AG':
			accs.append(i+1)

	return dons, accs

# using index starting at 0
def short_introns(dons, accs, minin):

	for d, a in zip(dons, accs):
		intron_length = a - d + 1
		if intron_length < minin:
			return True

	return False

# using index starting at 0
def short_exons(dons, accs, flank, minex, seq):

	# check 5' first exon
	fexlen = dons[0] - flank
	if fexlen < minex:
		return True

	# check 3' last exon
	lexbeg = accs[-1] + 1
	lexend = len(seq) - flank - 1 
	lexlen = lexend - lexbeg + 1
	if lexlen < minex:
		return True

	# check interior exons
	for i in range(len(dons)-1):
		if dons[i+1] - accs[i] - 1 < minex:
			return True

	return False	

def get_exons(dsites, asites, flank, seq):

	exons = []
	exons.append((flank, dsites[0]-1))
	for i in range(1, len(dsites)):
		exbeg = asites[i-1] + 1
		exend = dsites[i] - 1
		exons.append((exbeg, exend))
	exons.append((asites[-1]+1, len(seq)-flank-1))

	return exons	

def get_introns(dsites, asites):

	introns = []	
	for d, a in zip(dsites, asites):	
		introns.append((d, a))

	return introns

def apc(dons, accs, maxs, minin, minex, flank, seq):
	
	# use .copy() when appending dictionaries to lists
	apc_isoforms = []

	apc_isoform = {
		'seq': '',
		'beg': '',
		'end': '',
		'exons': [],
		'introns': [],
		'score': 0
	}	
	
	trials = 0
	short_introns_exons = 0
	nsites = min(len(dons), len(accs), maxs)
	for n in range(1, nsites+1):
		for dsites in combinations(dons, n):
			for asites in combinations(accs, n):
				trials += 1
				if short_introns(dsites, asites, minin): continue
				if short_exons(dsites, asites, flank, minex, seq): continue
				apc_isoform['seq'] = seq
				apc_isoform['beg'] = flank
				apc_isoform['end'] = len(seq) - flank - 1
				apc_isoform['exons'] = get_exons(dsites, asites, flank, seq)
				apc_isoform['introns'] = get_introns(dsites, asites)
				apc_isoforms.append(apc_isoform.copy())	
	return apc_isoforms, trials
