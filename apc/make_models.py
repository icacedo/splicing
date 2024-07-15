import argparse
import gzip
import isomod as im
import openturns as ot 
import csv

parser = argparse.ArgumentParser(
	description='Generates len, MM, and PWM models for apc')

parser.add_argument('--len_limit', type=int, metavar='<int>',
	required=False, help='size limit for length model')
parser.add_argument('--dntxt', type=str, metavar='<file>',
	required=False, help='input text file with donor site sequences')
parser.add_argument('--actxt', type=str, metavar='<file>',
	required=False, help='input text file with acceptor site sequences')
parser.add_argument('--outdir', type=str, metavar='<directory>',
	required=False, help='output directory name')
parser.add_argument('-mm', action='store_true', help='make only mm')
parser.add_argument('-len', action='store_true', help='make only len')

args = parser.parse_args()

def fdist_params(exinseqs, nbins=None, pre=None, size_limit=None):

	exinlens = im.get_exinbins(exinseqs, nbins=None, pre=None)[2]
	
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

def len_tsv_write(data, a, b, g, size_limit, fp, outdir=None):
	
	exinlen_yscores, exinlen_yvalues = aml.memoize_fdist(
		data, a, b, g, size_limit, pre=6)

	path = fp.split('/')
	root, ext = path[-1].split('.')
	if outdir:
		filename = outdir + root + '_len' + '.tsv'
	else:
		filename = root + '_len' + '.tsv'

	with open(filename, 'w', newline='') as tsvfile:
		writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
		writer.writerow(['% EVD params: '+'a: '+str(a)+' b: '+str(b)+' g '
				   +str(g)])
		writer.writerow(['% len '+root+' P', root+' log2(P/expect)'])
		for i in range(len(exinlen_yscores)):
			writer.writerow([exinlen_yvalues[i], exinlen_yscores[i]])
	tsvfile.close()

def mm_tsv_write(exins, fp, outdir=None):

	exinmm_scores, exinmm_probs, order = aml.make_mm(exins)
	
	path = fp.split('/')
	root, ext = path[-1].split('.')
	if outdir:
		filename = outdir + root + '_mm' + '.tsv'
	else:
		filename = root + '_mm' + '.tsv'

	with open(filename, 'w', newline='') as tsvfile:
		writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
		writer.writerow(['% mm '+root+' '+str(order+1)+'mer', 'mm '+root+' P',
			'mm '+root+' log2(P/0.25)'])
		for key in exinmm_scores:
			writer.writerow([key + 'A', exinmm_probs[key][0], 
				exinmm_scores[key][0]])
			writer.writerow([key + 'C', exinmm_probs[key][1], 
				exinmm_scores[key][1]])
			writer.writerow([key + 'G', exinmm_probs[key][2], 
				exinmm_scores[key][2]])
			writer.writerow([key + 'T', exinmm_probs[key][3], 
				exinmm_scores[key][3]])
	tsvfile.close()

if args.len_limit:
	len_limit = args.len_limit
else:
	len_limit = None

if args.extxt and args.mm:
	exons = aml.read_txt_seqs(args.extxt)
	mm_tsv_write(exons, args.extxt, args.outdir)
elif args.extxt and args.len:
	exons = aml.read_txt_seqs(args.extxt)
	data, a, b, g, limit = aml.fdist_params(exons, size_limit=len_limit) #***
	len_tsv_write(data, a, b, g, limit, args.extxt, args.outdir) #***
	#len_tsv_write(exons, args.extxt, args.outdir)	
elif args.extxt:
	exons = aml.read_txt_seqs(args.extxt)
	data, a, b, g, limit = aml.fdist_params(exons, size_limit=len_limit) #***
	len_tsv_write(data, a, b, g, limit, args.extxt, args.outdir) #***
	#len_tsv_write(exons, args.extxt, args.outdir)
	mm_tsv_write(exons, args.extxt, args.outdir)
							
if args.intxt and args.mm:
	introns = aml.read_txt_seqs(args.intxt)
	mm_tsv_write(introns, args.intxt, args.outdir)
elif args.intxt and args.len:
	introns = aml.read_txt_seqs(args.intxt)
	data, a, b, g, limit = aml.fdist_params(introns, size_limit=len_limit) #***
	len_tsv_write(introns, a, b, g, limit, args.intxt, args.outdir) #***
	#len_tsv_write(introns, args.intxt, args.outdir)
elif args.intxt:
	introns = aml.read_txt_seqs(args.intxt)
	data, a, b, g, limit = aml.fdist_params(introns, size_limit=len_limit) #***
	len_tsv_write(data, a, b, g, limit, args.intxt, args.outdir) #***
	#len_tsv_write(introns, args.intxt, args.outdir)
	mm_tsv_write(introns, args.intxt, args.outdir)

####################################################

def pdread(dictionary):

	for site in dictionary:
		A = None
		C = None
		G = None
		T = None
		for key in site:
			if key == 'A':
				A = site[key]
			if key == 'C':
				C = site[key]
			if key == 'G':
				G = site[key]
			if key == 'T':
				T = site[key]
		yield [A, C, G, T]

def de(num, pre=6):
	
	num2 = f'{num:.{pre}f}'
	return num2

def pwm_tsv_write(donacc, fp, outdir=None):

	pwm, ppm = aml.make_pwm(donacc)

	path = fp.split('/')
	root, ext = path[-1].split('.')
	if outdir:
		filename = outdir + root + '_pwm' + '.tsv'
	else:
		filename = root + '_pwm' + '.tsv'

	with open(filename, 'w', newline='') as tsvfile:
		writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
		writer.writerow(['% pwm '+root+' log2(PA/0.25)', 'log2(PC/0.25)', \
			'log2(PG/0.25)', 'log2(PT/0.25)'])
		for site in pdread(pwm):
			a = de(site[0])
			c = de(site[1])
			g = de(site[2])
			t = de(site[3])	
			writer.writerow([a, c, g, t])
		writer.writerow(['% ppm '+root+' PA', 'PC', 'PG', 'PT'])
		for site in pdread(ppm):
			a = de(site[0])
			c = de(site[1])
			g = de(site[2])
			t = de(site[3])
			writer.writerow([a, c, g, t])
		tsvfile.close()

if args.dntxt:
	donors = aml.read_txt_seqs(args.dntxt)
	pwm_tsv_write(donors, args.dntxt, args.outdir)

if args.actxt:
	acceptors = aml.read_txt_seqs(args.actxt)
	pwm_tsv_write(acceptors, args.actxt, args.outdir)







