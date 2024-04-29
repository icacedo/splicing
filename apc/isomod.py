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
	
	return new_ys
		
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