import sys
import gzip
import openturns as ot
import math

# unit testing in python
# https://www.dataquest.io/blog/unit-tests-python/

# for testing
fp = sys.argv[1]

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

print(memoize_fdist(fp))



########################################
##### End Length Model Section #########
########################################








