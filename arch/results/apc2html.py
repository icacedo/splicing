#!/usr/bin/env python3

import argparse
import glob
import os
import re
import sys
from operator import itemgetter

from grimoire.genome import Reader

def draw_transcript(tx, y):
	body = []
	if tx.score == '.': text = 'WormBase'
	else:               text = tx.score
	body.append(f'<text x="20" y="{y+10}" class="sm">{text}</text>')
	for exon in tx.exons:
		body.append(draw_exon(exon, y))
	for intron in tx.introns:
		body.append(draw_intron(intron, y))
	for cds in tx.cdss:
		body.append(draw_cds(cds, y))
	for utr in tx.utr5s:
		body.append(draw_utr(utr, y))
	for utr in tx.utr3s:
		body.append(draw_utr(utr, y))
	return '\n'.join(body)

def draw_exon(f, y):
	x = f'x="{f.beg}"';
	y = f'y="{y}"'
	w = f'width="{f.length}"'
	h = 'height="12"'
	if f.strand == '+': c = 'fill="#00f"'
	else:               c = 'fill="#f0f"'
	return f'<rect {x} {y} {w} {h} {c}><title>{f.beg}-{f.end}</title></rect>'

def draw_intron(f, y):
	x1 = f'x1="{f.beg}"';
	y1 = f'y1="{y+6}"'
	x2 = f'x2="{f.end}"'
	y2 = f'y2="{y+6}"'
	if f.strand == '+': c = 'stroke="#00f"'
	else:               c = 'stroke="#f0f"'
	return f'<line {x1} {y1} {x2} {y2} stroke-width="3" {c}><title>{f.seq_str()[0:5]}-{f.seq_str()[-6:]}</title></line>'

def draw_cds(f, y):
	x = f'x="{f.beg+1}"';
	y = f'y="{y+1}"'
	w = f'width="{f.length-1}"'
	h = 'height="10"'
	c = 'fill="#ff0"'
	return f'<rect {x} {y} {w} {h} {c}><title>{f.beg}-{f.end}</title></rect>'

def draw_utr(f, y):
	x = f'x="{f.beg+1}"';
	y = f'y="{y+1}"'
	w = f'width="{f.length-1}"'
	h = 'height="10"'
	c = 'fill="brown"'
	return f'<rect {x} {y} {w} {h} {c}><title>{f.beg}-{f.end}</title></rect>'

def draw_splice(f, y):
	x = f'x="{f.beg}"';
	y = f'y="{y}"'
	w = f'width="{f.length}"'
	h = 'height="5"'
	if f.strand == '+': c = 'fill="#00f"'
	else:               c = 'fill="#f0f"'
	return f'<rect {x} {y} {w} {h} {c}><title>{f.score}</title></rect>'

def get_isoforms(fasta,  params):
	os.system(f'isoformer {fasta} {params} > tmp.gff')
	genome = Reader(gff='tmp.gff', fasta=fasta)
	chrom = next(genome)
	for gene in chrom.ftable.build_genes():
		for tx in gene.transcripts():
			tx.use_longest_orf()
			yield tx

# Hard-coded badness
apc = 'apc'      # dir of fasta and gff
tsv = 'apc.tsv'  # file of weights
out = 'html'     # outputdir
mod = '../data'  # directory where models are

# Isoformer params
params = {}
fitness = {}
with open(tsv) as fp:
	header = fp.readline()
	for line in fp.readlines():
		seq, don, acc, emm, imm, elen, ilen, cost, fit = line.split()
		params[seq] = [
			'--dpwm',  f'{mod}/donor.pwm',    f'--wdpwm {don}',
			'--apwm',  f'{mod}/acceptor.pwm', f'--wapwm {acc}',
			'--emm',   f'{mod}/exon.mm',      f'--wemm {emm}',
			'--imm',   f'{mod}/intron.mm',    f'--wimm {imm}',
			'--elen',  f'{mod}/exon.len',     f'--welen {elen}',
			'--ilen',  f'{mod}/intron.len',   f'--wilen {ilen}',
			'--icost', f'{cost}'
		]
		fitness[seq] = fit

# Text
style = '''
<style>
  .sm { font: 12px sans-serif; }
  .lg { font: 15px sans-serif; }
</style>
'''

# Main
for ff in glob.glob(f'{apc}/*.fa'):
	gf = ff[:-2] + 'gff3'
	ch = ff[4:-3]
	print(f'building html for {ch}', file=sys.stderr)

	genome = Reader(gff=gf, fasta=ff)
	chrom = next(genome)

	body = [] # body of SVG
	
	## Title
	off = 20 # offset
	body.append(f'<rect fill="#ddd" x="0" y="0" width="100%" height="100%"/>')
	body.append(f'<text x="20" y="20" class="lg">{chrom.name} {chrom.desc}</text>')

	## Transcripts
	for gene in chrom.ftable.build_genes():
		for tx in gene.transcripts():
			off += 20
			body.append(draw_transcript(tx, off))

	## Fitmess
	body.append(f'<text x="20" y="{off+30}" class="sm">Fitness: {fitness[ch]}</text>')

	## RNASeq_splice (should scale by score)
	rss = []
	for f in chrom.ftable.features:
		if f.source == 'RNASeq_splice':
			rss.append(f)
	off += 20
	for splice in sorted(rss, key=lambda x: x.score, reverse=True):
		off += 5
		body.append(draw_splice(splice, off))

	## Predicted isoforms & CDS
	for iso in get_isoforms(ff, ' '.join(params[ch])):
		off += 20
		body.append(draw_transcript(iso, off))
	

	xmlns = 'xmlns="http://www.w3.org/2000/svg"'
	vbox = f'viewBox="0 0 {len(chrom.seq)} {off+40}"'
	
	## html goodies
	match = re.search('(WBGene\d+)', chrom.desc)
	wbgene = match.group(1)
	
	with open(f'{out}/{chrom.name}.html', 'w') as fp:
		fp.write(f'<html><title>{chrom.name}</title><body>\n')
		fp.write(f'<h1>{chrom.name} {len(chrom.seq)} bp</h1>\n')
		fp.write(f'<a href="https://wormbase.org/species/c_elegans/gene/{wbgene}">{wbgene}</a>\n')
		fp.write(f'<svg {vbox} {xmlns}>\n')
		fp.write(style)

		for line in body:
			fp.write(line)
			fp.write('\n')
		fp.write('</svg>\n')
		
		fp.write('</body></html>')

	#sys.exit()