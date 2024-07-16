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
parser.add_argument('--len_limit', type=int, default=1000, metavar='<int>',
	required=False, help='size limit for length model %(default)d')
parser.add_argument('--outdir', type=str, metavar='<directory>',
	required=False, help='output directory name')

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

def write_len(data, a, b, g, size_limit, fp, outdir=None):
	
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

#eseqs, iseqs, dseqs, aseqs = im.get_all_tn_seqs(args.wb_dir)

def get_gff_tn_seqs(seq, gff):

	eseqs = []
	with open(gff, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			if line[1] == 'WormBase':
				print(line)

gffs = {}
fastas = {}
for file in os.listdir(args.wb_dir):
	id = file.split('.')[1]
	if file.endswith('gff3'):
		gffs[id] = f'{args.wb_dir}{file}'
	if file.endswith('fa'):
		fastas[id] = f'{args.wb_dir}{file}'

seq = im.read_fasta(fastas['3068'])
print(seq[0])

get_gff_tn_seqs(seq[1], gffs['3068'])

'''
for gene in gffs:
	seq = im.read_fasta(fastas[gene])
	#eseqs, iseqs, dseqs, aseqs = get_gff_tn_seqs(seq[1], gffs[gene])
	eseqs = get_gff_tn_seqs(seq[1], gffs[gene])
	print(eseqs)
	#print(iseqs)
	#print(dseqs)
	#print(aseqs)
	print(seq[0])
	break
'''

# testing isoforms/modelbuilder

import genome


exons = []
genome = genome.Reader(gff=gffs['3068'], fasta=fastas['3068'])
tx = next(genome).ftable.build_genes()[0].transcripts()[0]
for f in tx.exons: exons.append(f.seq_str())

print(exons)

e, i, d, a = im.get_gff_tn_seqs(seq[1], gffs['3068'])
print(e)

# ch.6441 overlaps with another gene
# ch.3068 is on the + strand on WB
# ch.862 is on the - strand on WB
# should the models be trained on only the canonical sequences?
# why does the exon include the 5' utr?

'''
accs = []
dons = []
exons = []
introns = []
for ff in glob.glob(f'{args.wb_dir}/*.fa'):
	gf = ff[:-2] + 'gff3'
	genome = genome.Reader(gff=gf, fasta=ff)
	tx = next(genome).ftable.build_genes()[0].transcripts()[0]
	for f in tx.exons: exons.append(f.seq_str())
	for f in tx.introns:
		iseq = f.seq_str()
		dons.append(iseq[0:5])
		accs.append(iseq[-6:])
		introns.append(iseq)
'''