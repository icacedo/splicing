import argparse
import os
import gzip
import isomod as im
import openturns as ot 
import csv
import glob

parser = argparse.ArgumentParser(description='Generates len, MM, and PWM \
	models for apc based on sequences in the apc dataset')
parser.add_argument('wb_dir', type=str, metavar='<directory>', 
	help='directory with apc dataset gff and fasta files')
parser.add_argument('--outdir', type=str, metavar='<directory>',
	required=False, help='output directory name')

args = parser.parse_args()

def fdist_params(exinseqs, len_limit):

	exinlens = im.get_exinbins(exinseqs)[0]
	
	if len_limit:
		sample = ot.Sample([[x] for x in exinlens if x < len_limit])
	else:
		len_limit = max(exinlens)
		sample = ot.Sample([[x] for x in exinlens if x < len_limit])	

	distFrechet = ot.FrechetFactory().buildAsFrechet(sample)

	a = distFrechet.getAlpha()
	b = distFrechet.getBeta()
	g = distFrechet.getGamma()
	
	return exinlens, a, b, g

gffs = {}
fastas = {}
for file in os.listdir(args.wb_dir):
	id = file.split('.')[1]
	if file.endswith('gff3'):
		gffs[id] = f'{args.wb_dir}{file}'
	if file.endswith('fa'):
		fastas[id] = f'{args.wb_dir}{file}'

exons = []
introns = []
dons = []
accs = []
for gid in gffs:
	seq = im.read_fasta(fastas[gid])
	subseqs = im.get_subseqs(seq[1], gffs[gid])
	for e in subseqs[0]: exons.append(e)
	for i in subseqs[1]: introns.append(i)
	for d in subseqs[2]: dons.append(d)
	for a in subseqs[3]: accs.append(a)

elens, a, b, g = fdist_params(exons, 1000)

v = im.memoize_fdist(elens, a, b, g, 25, 1000)

print(v)


import math

def create_length_model(seqs, lmin, lmax):
	# get the lengths
	lens = [len(seq) for seq in seqs]

	# train Frechet
	sample = ot.Sample([[x] for x in lens if x < lmax])
	f = ot.FrechetFactory().buildAsFrechet(sample)
	a = f.getAlpha()
	b = f.getBeta()
	g = f.getGamma()
	
	# create histogram from frechet
	pdf = []
	for x in range(lmax):
		if x < g: pdf.append(0)
		else:
			z = (x-g)/b
			pdf.append((a/b) * z**(-1-a) * math.exp(-z**-a))

	# create leading zeros and rescale
	for i in range(lmin): pdf[i] = 0
	total = sum(pdf)
	for i in range(len(pdf)): pdf[i] /= total	

	return pdf

pdf = create_length_model(exons, 25, 1000)
pdf2 = []
for p in pdf:
	pdf2.append(f'{p:.{6}f}')

print(pdf2)

if v == pdf2: print('wow')







# testing code from isoforms/modelbuilder
'''
import genome

exons = []
introns = []
dons = []
accs = []
genome = genome.Reader(gff=gffs['3068'], fasta=fastas['3068'])
tx = next(genome).ftable.build_genes()[0].transcripts()[0]
for f in tx.exons: exons.append(f.seq_str())
for f in tx.introns:
	iseq = f.seq_str()
	dons.append(iseq[0:5])
	accs.append(iseq[-6:])
	introns.append(iseq)

s = im.read_fasta(fastas['3068'])
subs = im.get_subseqs(s[1], gffs['3068'])

if exons == subs[0]: print('same e')
if introns == subs[1]: print('same i')
if dons == subs[2]: print('same d')
if accs == subs[3]: print('same a')


'''
# ch.6441 overlaps with another gene
# ch.3068 is on the + strand on WB
# ch.862 is on the - strand on WB


