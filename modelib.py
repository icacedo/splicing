import sys
import gzip
import openturns as ot
import math
from itertools import combinations

# unit testing in python
# https://www.dataquest.io/blog/unit-tests-python/

# for testing
#fp = sys.argv[1]

########################################
##### Begin File Reading Section #######
########################################

def read_fastas(fastafile):

	with open(fastafile) as fp:

		seqid = ''
		seq = ''	

		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('>'): seqid += line
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
def read_txt_seqs(sfile):
	
	exin_seqs = []
	with open(sfile, 'r') as fp:	
		exin_seq = ''
		for line in fp.readlines():
			line = line.rstrip()
			exin_seq = line
			exin_seqs.append(exin_seq)
	return exin_seqs
	fp.close()

########################################
##### End File Reading Section #########
########################################

########################################
##### Begin Length Model Section #######
########################################

def read_intfile(fp):

	if fp.endswith(".gz"):
		with gzip.open(fp, 'r') as intfile:            
			for line in fp.readlines():
				line = line.rstrip()
				if isinstance(line, bytes):
					line = line.decode()
				yield line

	else:
		with open(fp, 'r') as intfile:
			for line in intfile.readlines():
				line = line.rstrip()
				yield line

# returns a count/frequency  histogram or raw counts
def get_intbins(fp, nbins=None, prec=None):
	
	intlines = []
	total_obs = 0
	for l in read_intfile(fp):
		intlines.append(l)
		total_obs += 1

	intsizes = []
	for intron in intlines:
		intsizes.append(len(intron))
	
	intfreqs = []
	total = len(intsizes)
	for count in intsizes:
		fq = count/total
		if prec:
			fq2 = f"{fq:.{prec}f}"
			intfreqs.append(float(fq2))
		else: intfreqs.append(float(fq))
	
	intcount_bins = [0 for x in range(max(intsizes)+1)]
	for i in range(len(intsizes)):
		intcount_bins[intsizes[i]] += 1

	intfreq_bins = []
	for i in range(len(intcount_bins)):
		fqb = intcount_bins[i]/total_obs
		if nbins:
			fqb2 = f"{fqb:.{prec}f}"
			intfreq_bins.append(float(fqb))
		else: intfreq_bins.append(float(fqb))	

	# only append up to nbins, for testing
	return intcount_bins[:nbins], intfreq_bins[:nbins], intsizes, intfreqs

##### smoothing ########################

# definitions of rectangular and triangular smoothing found here:
# https://terpconnect.umd.edu/~toh/spectrum/Smoothing.html

# rectangular smoothing
def rec_smoo(intbins, m=5, prec=None):
	
	smoodata = []
	for i in range(len(intbins)):
		m2 = int((m/2) + 0.5 - 1)
		bef = intbins[i-m2:i]
		now = intbins[i]
		aft = intbins[i+1:i+m2+1]
		total = sum(bef) + now + sum(aft)
		smoopt = total/m
		if prec:
			fmat = f'{smoopt:.{prec}f}'
			smoodata.append(fmat)
		else:
			smoodata.append(smoopt)
		#print(bef, now, aft, total, smoopt)
	
	return smoodata

# triangular smoothing
def tri_smoo(intbins, m=5, prec=None):	

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
		if prec: 
			fmat = f'{smoopt:.{prec}f}'
			smoodata.append(fmat)
		else:
			smoodata.append(smoopt)
		#print(bef, now, aft, total, smoopt)

	return smoodata

##### curve fitting ####################

def frechet_pdf(x, a, b, g):
	if x < g: return 0
	z = (x-g)/b
	term1 = (a/b)
	term2 = z**(-1-a)
	term3 = math.exp(-z**-a)
	return term1 * term2 * term3

def memoize_fdist(fp, nbins=None, prec=None, size_limit=250):

	data = get_intbins(fp, nbins=None, prec=None)[2]
	
	# data is used to get the parameters
	# this function does not score the data
	# only add introns below certain size to sample 
	sample = ot.Sample([[x] for x in data if x < size_limit])

	distFrechet = ot.FrechetFactory().buildAsFrechet(sample)

	a = distFrechet.getAlpha()
	b = distFrechet.getBeta()
	g = distFrechet.getGamma()
		
	x_values = []
	y_values = []
	y_scores = []
	for i in range(min(len(data), size_limit)):
		x_values.append(i)
		y = frechet_pdf(i, a, b, g)
		y_values.append(y)
		if y == 0: y_scores.append(-99)
		else: y_scores.append(math.log2(y))

	data = {
		'x': x_values,
		'y': y_values
	}
	
	# only scores are useful?
	return y_scores

#print(memoize_fdist(fp))



########################################
##### End Length Model Section #########
########################################

########################################
##### Begin Markov Model Section #######
########################################

def make_mm(exin_seqs, order=3):

	order = 3
	context = {}
	for seq in exin_seqs:
		for i in range(len(seq)-order):
			prev = seq[i:i+order]
			now = seq[i+order]
			if prev not in context:
				context[prev] = now
			else:
				context[seq[i:i+order]] += now
	mm = {}
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
		if nts not in mm:
			mm[nts] = (A/d, C/d, G/d, T/d)
	return mm

########################################
##### End Markov Model Section #########
########################################

########################################
##### Begin APC Section ################
########################################

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
				print(apc_isoform)
	print(trials)








