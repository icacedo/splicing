import sys
import gzip
import math
import os
import openturns as ot
from itertools import combinations

################################
##### File Reading Section #####
################################

def read_fasta(fastafile):

	with open(fastafile) as fp:

		seqid = ''
		seq = ''	

		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('>'): seqid = line[1:]
			else: seq += line
		return (seqid, seq)			
	fp.close()	

def read_gff_sites(seq, gff, gtag=True):

	dons = []
	accs = []

	with open(gff) as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'intron':
				beg = int(line[3]) - 1
				end = int(line[4])
				if gtag:
					if seq[beg:beg+2] == 'GT': 
						dons.append(beg-1)
					if seq[end-2:end] == 'AG': 
						accs.append(end-1)
				if not gtag: 
					dons.append(beg-1)
					accs.append(end-1)
	fp.close()	

	return sorted(set(dons)), sorted(set(accs))

# get exon/intron/donor/acceptor training seqs from gffs
'''
def get_gff_tn_seqs(seq, gff):

	iseqs = []
	eseqs = []
	dseqs = []
	aseqs = []
	with open(gff, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			if line[2] == 'intron': 
				beg = int(line[3])
				end = int(line[4])
				iseq = seq[beg-1:end]
				if len(iseq) <= 20:
					print(gff, iseq)
				iseqs.append(iseq)
				dseq = seq[beg-1:beg+4]
				dseqs.append(dseq)
				aseq = seq[end-6:end]
				aseqs.append(aseq)
			if line[2] == 'exon':
				beg = int(line[3])
				end = int(line[4])
				eseq = seq[beg-1:end]
				if len(eseq) <= 20:
					print(gff, eseq)
				eseqs.append(eseq)

	return eseqs, iseqs, dseqs, aseqs
'''
# work in progress
def get_top_exins(seq, gff):

	wb_ints = []
	wb_exos = []
	rsplice_ints = {}
	with open(gff, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip().split('\t')
			if line[1] == 'WormBase':
				if line[2] == 'intron':
					intron = line[3], line[4]
					wb_ints.append(intron)
				if line[2] == 'exon':
					exon = line[3], line[4]
					wb_exos.append(exon)
			if line[1] == 'RNASeq_splice':
				if line[2] == 'intron':
					intron = int(line[3]), int(line[4])
					rsplice_ints[intron] = float(line[5])

	rsplice_ints = {i: j for i, j in sorted(rsplice_ints.items(), 
								 key=lambda tem: tem[1], reverse=True)}
	total = 0
	for intron in rsplice_ints:
		total += rsplice_ints[intron]

	top_ints = []
	scores = 0
	for intron in rsplice_ints:
		if scores >= 99: break
		top_ints.append(intron)
		scores += (rsplice_ints[intron]/total)*100
	
	print(sorted(top_ints), '#$#')
	print(wb_exos)
	exends = []
	exbegs = []
	for intron in sorted(top_ints):
		exends.append(intron[0]-1)
		exbegs.append(intron[1]+1)
	exends = sorted(exends)
	exbegs = sorted(exbegs)
	print(exends)
	print(exbegs)
	for i in range(len(exends)-1):
		print(i, exbegs[i], exends[i+1])
		
def get_all_tn_seqs(wb_dir):

	fastas = {}
	gffs = {}
	for fname in os.listdir(wb_dir):
		fid = fname.split('.')[1]
		if fname.endswith('.fa'):
			fastas[fid] = wb_dir + fname
		if fname.endswith('.gff3'):
			gffs[fid] = wb_dir + fname
			
	pairs = {}
	for fid in fastas:
		pairs[fid] = [fastas[fid], gffs[fid]]

	all_eseqs = []
	all_iseqs = []
	all_dseqs = []
	all_aseqs = []
	for fid in pairs:
		seqid, seq = read_fasta(pairs[fid][0])
		gff = pairs[fid][1]
		eseqs, iseqs, dseqs, aseqs = get_gff_tn_seqs(seq, gff)
		for e, i, d, a in zip(eseqs, iseqs, dseqs, aseqs):
			all_eseqs.append(e)
			all_iseqs.append(i)
			all_dseqs.append(d)
			all_aseqs.append(a)

	return all_eseqs, all_iseqs, all_dseqs, all_aseqs

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

def frechet_pdf(x, a, b, g):

	if x < g: return 0
	z = (x-g)/b
	term1 = (a/b)
	term2 = z**(-1-a)
	term3 = math.exp(-z**-a)
	return term1 * term2 * term3

def fdist_params(exinseqs, nbins=None, pre=None, size_limit=None):

	exinlens = get_exinbins(exinseqs, nbins=None, pre=None)[2]
	
	if size_limit:
		sample = ot.Sample([[x] for x in exinlens if x < size_limit])
	else:
		size_limit = max(exinlens)
		sample = ot.Sample([[x] for x in exinlens if x < size_limit])	

	distFrechet = ot.FrechetFactory().buildAsFrechet(sample)

	a = distFrechet.getAlpha()
	b = distFrechet.getBeta()
	g = distFrechet.getGamma()
	
	return exinlens, a, b, g, size_limit

# at size_limit 500-max exon/intron fit frechet dist
# lower than 500, no longer frechet dist
# pre stands for precision (decimals)
def memoize_fdist(data, a, b, g, size_limit, pre=None):
		
	assert size_limit >= 500, "limit too small for frechet"

	x_values = []
	y_values = []
	y_scores = []
	# is expectation necessary for extreme value distribution?
	expect = 1/size_limit
	for i in range(min(len(data), size_limit)):
		x_values.append(i)
		y = frechet_pdf(i, a, b, g)
		if pre:
			y2 = f"{y:.{pre}f}"
			y_values.append(y2)
			if y == 0: y_scores.append(-100)
			else:
				ys = math.log2(y/expect)
				ys2 =  f"{ys:.{pre}f}"
				y_scores.append(ys2)
		else: 
			y_values.append(y)
			if y == 0: y_scores.append(-100)
			else: y_scores.append(math.log2(y/expect))			
	
	data = {
		'x': x_values,
		'y': y_values
	}
	
	return y_scores, y_values

# assign x values (exon/intron lengths) below smallest observed 0 prob
# use this to make the model files
def zero_prob(exinlens, y_values, pre=6):

	small = min(exinlens)
	y_vals = []
	total = 0
	for i in range(len(y_values)):
		if i < small:
			y_vals.append(0)
		else:
			total += y_values[i]
			y_vals.append(y_values[i])
	new_ys = [float(f'{x/total:.{pre}f}') for x in y_vals]
	'''
if args.elen:
	re_elen = im.read_len(args.elen)
else:
	re_elen = None
'''
	return new_ys

###########################
##### Scoring Section #####	
###########################

# len model scoring

def read_len(len_model):

	re_len = []
	with open(len_model, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('%'): continue
			re_len.append(float(line))

	return re_len

def score_len(re_len, exin):

	if re_len == None: return 0

	length = exin[1] - exin[0]
	len_prob = re_len[length]
	if len_prob == 0:
		len_score = -99
	else:
		exp = 1/len(re_len)
		len_score = math.log2(len_prob/exp)

	return len_score


# mm scoring

def score_mm(re_mm, exin, seq, dpwm=None, apwm=None):

	beg = exin[0]
	end = exin[1] + 1
	exin_seq = seq[beg:end]

	k = 0
	for kmer in re_mm:
		k = len(kmer)
		break
	
	if dpwm and apwm:
		exin_seq = exin_seq[len(dpmw):-len(apwm)]

	mm_score = 0
	for i in range(len(exin_seq)):
		if len(exin_seq[i:i+k]) == k:
			kmer = exin_seq[i:i+k]
			mm_score += math.log2(float(re_mm[kmer])/0.25)

	return mm_score

# pwm model scoring

def get_daseq(intron, seq):

	d_start = intron[0]
	d_end = d_start + 5
	a_end = intron[1] + 1
	a_start = a_end -6
	d_seq = seq[d_start:d_end]
	a_seq = seq[a_start:a_end]
	
	return d_seq, a_seq

def score_pwm(daseq, pwm):
	
	da_score = 0
	count = 0
	for i in range(len(daseq)):
		if daseq[i] == 'A':
			da_score += math.log2(float(pwm[count][0])/0.25)
		if daseq[i] == 'C':
			da_score += math.log2(float(pwm[count][1])/0.25)
		if daseq[i] == 'G':
			da_score += math.log2(float(pwm[count][2])/0.25)
		if daseq[i] == 'T':
			da_score += math.log2(float(pwm[count][3])/0.25)
		count += 1

	return da_score

def get_entropy(probs):

	h = 0
	for p in probs:
		h -= p * math.log2(p)
	return h













##### Markov Model scoring #####

def read_mm(mm_model):

	mm_scores = {}
	with open(mm_model, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('%'): continue
			if line == '': continue
			line = line.split(' ')
			mm_scores[line[0]] = float(line[1])

	return mm_scores

##### PWM model scoring #####

def read_pwm(pwm_model):

	pwm_scores = []
	with open(pwm_model, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('%'): continue
			line = line.split(' ')
			pwm_scores.append(line)
	
	return pwm_scores








#######################################
##### All Biological Combinations #####
#######################################

def get_gtag(seq, flank, minex):

	dons = []
	accs = []
	for i in range(flank + minex, len(seq) - flank - minex - 1):
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

def abc(dons, accs, maxs, minin, minex, flank, seq):
	
	# use .copy() when appending dictionaries to lists
	abc_isoforms = []

	abc_isoform = {
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
				abc_isoform['seq'] = seq
				abc_isoform['beg'] = flank
				abc_isoform['end'] = len(seq) - flank - 1
				abc_isoform['exons'] = get_exons(dsites, asites, flank, seq)
				abc_isoform['introns'] = get_introns(dsites, asites)
				abc_isoforms.append(abc_isoform.copy())	
	return abc_isoforms, trials






