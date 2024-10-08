# get all apc isoforms once, re-run the scoring part many times to get icost
import argparse
import pickle
import os
import sys
import apc_model_lib as aml

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
	help='input single sequence fasta file')
parser.add_argument('--gff', type=str, metavar='<file>', required=False,
	help='input .gff3 for single gene')
parser.add_argument('--outdir', type=str, metavar='<outdir path>',
	required=True, help='/path/ to outdir')

# apc parameters
parser.add_argument('--max_splice', required=False, type=int, default=3,
	metavar='<int>', help='maximum number of splicing events %(default)d')
parser.add_argument('--min_intron', required=False, type=int, default=25,
	metavar='<int>', help='minimum length of intron %(default)d')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length of exon %(default)d')
parser.add_argument('--flank', required=False, type=int, default=100,
	metavar='<int>', help='length of genomic flank on each side %(default)d')

args = parser.parse_args()

seqid = None
seq = None
for seqid, seq in aml.read_fastas(args.fasta):
	seqid = seqid
	seq = seq

if args.gff:
	dons, accs = aml.read_gff_sites(seq, args.gff)
else:
	dons, accs = aml.get_gtag(seq)

maxs = args.max_splice
minin = args.min_intron
minex = args.min_exon
flank = args.flank

apc_isoforms, trials = aml.apc(dons, accs, maxs, minin, minex, flank, seq)

fpath = args.fasta
fname = fpath.split('/')[-1]
ID = fname.split('.')[1]
outdir = args.outdir
os.makedirs(os.path.dirname(outdir), exist_ok=True)

if args.gff:
	name = outdir+'ch.'+ID+'.apc_isoforms.bli.pkl'
else:
	name = outdir+'ch.'+ID+'.apc_isoforms.pkl'
print(name)
with open(name, 'wb') as pick:
	pickle.dump(apc_isoforms, pick)

with open(name, 'rb') as pick:
	pickell = pickle.load(pick)

assert apc_isoforms == pickell, 'pickled incorrectly'
