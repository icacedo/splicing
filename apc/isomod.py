import sys
import gzip
import math
import os
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

	return sorted(set(dons)), sorted(set(accs))

def get_subseqs(seq, gff):

	eseqs = []
	iseqs = []
	dseqs = []
	aseqs = []
	with open(gff, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			if line[1] == 'WormBase' and line[2] == 'exon':
				beg = int(line[3])
				end = int(line[4])
				eseq = seq[beg-1:end]
				eseqs.append(eseq)
			if line[1] == 'WormBase' and line[2] == 'intron':
				beg = int(line[3])
				end = int(line[4])
				iseq = seq[beg-1:end]
				iseqs.append(iseq)
				dseq = seq[beg-1:beg+4]
				dseqs.append(dseq)
				aseq = seq[end-6:end]
				aseqs.append(aseq)

	return [eseqs, iseqs, dseqs, aseqs]

# work in progress
# wormbase gffs have more than the canonical exon/intron/donor/acceptor
# all sequences are considered equal when training feature models
# would it be worth it to only train on the features with the most transcripts?
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

# returns a count/frequency histogram or raw counts
# uses a list of exons/introns as input
def get_exinbins(exinseqs):

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
		fq2 = f"{fq:.{6}f}"
		freqs.append(float(fq2))
	
	count_bins = [0 for x in range(max(sizes)+1)]
	for i in range(len(sizes)):
		count_bins[sizes[i]] += 1

	freq_bins = []
	for i in range(len(count_bins)):
		fqb = count_bins[i]/total_obs	
		fqb2 = f"{fqb:.{6}f}"
		freq_bins.append(float(fqb2))	

	return sizes, freqs

def frechet_pdf(x, a, b, g):

	if x < g: return 0
	z = (x-g)/b
	term1 = (a/b)
	term2 = z**(-1-a)
	term3 = math.exp(-z**-a)
	return term1 * term2 * term3

# at size_limit 500-max exon/intron fit frechet dist
# lower than 500, no longer frechet dist
def memoize_fdist(data, a, b, g, minlen, maxlen):

	assert maxlen >= 500, "max length too small for frechet"

	yvals = []
	exp = 1/maxlen
	for x in range(min(len(data), maxlen)):
		y = frechet_pdf(x, a, b, g)
		yvals.append(y)

	# assign x values below size limit to 0
	yvals2 = []
	total = 0
	for i in range(len(yvals)):
		if i < minlen:
			yvals2.append(0)
		else:
			yvals2.append(yvals[i])
			total += yvals[i]

	for i in range(len(yvals2)):
		yvals2[i] = f"{yvals2[i]/total:.{6}f}"

	return yvals2

def len_write(data, fname, outdir=None):

	assert fname == 'exon' or fname == 'intron', 'file name not valid'

	if outdir:
		fpath = outdir + f"{fname}.len"
	else:
		fpath = f"{fname}.len"

	with open(fpath, 'w', newline='') as file:
		file.write(f"% LEN models/{fname}.len {len(data)}\n")
		for prob in data:
			file.write(f"{prob}\n")

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
				context[prev] += now

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
		mm_probs[nts] = [float(f"{x:.6f}") for x in [A/d, C/d, G/d, T/d]]		

	return mm_probs

def mm_write(data, fname, outdir=None):

	assert fname == 'exon' or fname == 'intron', 'file name not valid'

	if outdir:
		fpath = outdir + f"{fname}.mm"
	else:
		fpath = f"{fname}.mm"

	with open(fpath, 'w', newline='') as file:
		file.write(f"% MM models/{fname}.mm {len(data)*4}\n")
		alph = ['A', 'C', 'G', 'T']
		for cxt in data:
			for i in range(len(data[cxt])):
				file.write(f"{cxt}{alph[i]} {data[cxt][i]}\n")
			if cxt != list(data)[-1]:
				file.write('\n')

################################
##### PWM Model Section ########
################################

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

	for i in range(len(ppm)):
		for n in ppm[i]:
			ppm[i][n] = f"{ppm[i][n]:.6f}"	

	return ppm

def pwm_write(data, fname, outdir=None):

	assert fname == 'donor' or fname == 'acceptor', 'file name not valid'

	if outdir:
		fpath = outdir + f"{fname}.pwm"
	else:
		fpath = f"{fname}.pwm"

	with open(fpath, 'w', newline='') as file:
		file.write(f"% PWM models/{fname}.pwm {len(data)}\n")
		for site in data:
			file.write(f"{site['A']} {site['C']} {site['G']} {site['T']}\n")

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
	if length > len(re_len):
		len_prob == 0
	else:
		len_prob = re_len[length]
	if len_prob == 0:
		len_score = -100
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
		exin_seq = exin_seq[len(dpwm):-len(apwm)]

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
			if float(pwm[count][0]) == 0: da_score += -100
			else:
				da_score += math.log2(float(pwm[count][0])/0.25)
		if daseq[i] == 'C':
			if float(pwm[count][1]) == 0: da_score += -100
			else:
				da_score += math.log2(float(pwm[count][1])/0.25)
		if daseq[i] == 'G':
			if float(pwm[count][2]) == 0: da_score += -100
			else:
				da_score += math.log2(float(pwm[count][2])/0.25)
		if daseq[i] == 'T':
			if float(pwm[count][3]) == 0: da_score += -100
			else:
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
##### All Possible Combinations #####
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






