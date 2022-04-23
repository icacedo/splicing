#!/usr/bin/env python3

import argparse
import glob
import os
import sys
from operator import itemgetter

from grimoire.genome import Reader

def draw_transcript(tx, y):
	body = []
	for exon in tx.exons:
		body.append(draw_exon(exon, y))
	for intron in tx.introns:
		body.append(draw_intron(intron, y))
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
	return f'<rect {x} {y} {w} {h} {c}/>'

def draw_intron(f, y):
	x1 = f'x1="{f.beg}"';
	y1 = f'y1="{y+6}"'
	x2 = f'x2="{f.end}"'
	y2 = f'y2="{y+6}"'
	if f.strand == '+': c = 'stroke="#00f"'
	else:               c = 'stroke="#f0f"'
	return f'<line {x1} {y1} {x2} {y2} {c}/>'

def draw_utr(f, y):
	x = f'x="{f.beg}"';
	y = f'y="{y}"'
	w = f'width="{f.length}"'
	h = 'height="12"'
	c = 'fill="brown"'
	return f'<rect {x} {y} {w} {h} {c}><title>{f.score}</title></rect>'

def draw_splice(f, y):
	x = f'x="{f.beg}"';
	y = f'y="{y}"'
	w = f'width="{f.length}"'
	h = 'height="5"'
	if f.strand == '+': c = 'fill="#00f"'
	else:               c = 'fill="#f0f"'
	return f'<rect {x} {y} {w} {h} {c}><title>{f.score}</title></rect>'


# CLI
parser = argparse.ArgumentParser(
	description='apc image maker')
parser.add_argument('apc', type=str, metavar='<apc dir>',
	help='directory with apc set')
parser.add_argument('out', type=str, metavar='<out dir>',
	help='output directory')
arg = parser.parse_args()

# Text
style = '''
<style>
  .sm { font: 12px sans-serif; }
  .lg { font: 15px sans-serif; }
</style>
'''

for ff in glob.glob(f'{arg.apc}/*.fa'):
	gf = ff[:-2] + 'gff3'

	genome = Reader(gff=gf, fasta=ff)
	chrom = next(genome)

	body = [] # body of SVG
	off = 20 # offset
	body.append(f'<rect fill="#ddd" x="0" y="0" width="100%" height="100%"/>')
	body.append(f'<text x="20" y="20" class="lg">{chrom.name} {chrom.desc}</text>')

	for gene in chrom.ftable.build_genes():
		for tx in gene.transcripts():
			off += 20
			body.append(draw_transcript(tx, off))

	rss = []
	for f in chrom.ftable.features:
		if f.source == 'RNASeq_splice':
			rss.append(f)
	off += 20
	for splice in sorted(rss, key=lambda x: x.score, reverse=True):
		off += 5
		body.append(draw_splice(splice, off))

	xmlns = 'xmlns="http://www.w3.org/2000/svg"'
	vbox = f'viewBox="0 0 {len(chrom.seq)} {500}"'
	
	with open(f'{arg.out}/{chrom.name}.svg', 'w') as fp:
		fp.write(f'<svg {vbox} {xmlns}>\n')
		fp.write(style)

		for line in body:
			fp.write(line)
			fp.write('\n')
		fp.write('</svg>\n')

	#sys.exit()


	