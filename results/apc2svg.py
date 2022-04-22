#!/usr/bin/env python3

import argparse
import os
import sys

from grimoire.genome import Reader


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

# CLI
parser = argparse.ArgumentParser(
	description='apc image maker')
parser.add_argument('apc', type=str, metavar='<apc dir>',
	help='directory with apc set')
arg = parser.parse_args()


for ff in glob.glob(f'{arg.apc}/*.fa'):
	gf = ff[:-2] + 'gff3'

	genome = Reader(gff=gf, fasta=ff)
	chrom = next(genome)

"""

body = [] # body of SVG
off = 20 # offset



chrom = next(Reader(arg.fasta, arg.gff))



body.append(f'<rect fill="#ddd" x="0" y="0" width="100%" height="100%"/>')
body.append(f'<text x="20" y="20" class="lg">{chrom.name} {chrom.desc}</text>')

for gene in chrom.ftable.build_genes():
	for tx in gene.transcripts():
		off += 20
		for exon in tx.exons:
			body.append(draw_exon(exon, off))
		for intron in tx.introns:
			body.append(draw_intron(intron, off))


# SVG output
style = '''
<style>
  .sm { font: 12px sans-serif; }
  .lg { font: 15px sans-serif; }
</style>
'''

xmlns = 'xmlns="http://www.w3.org/2000/svg"'
vbox = f'viewBox="0 0 {len(chrom.seq)} {500}"'
print(f'<svg {vbox} {xmlns}>')
print(style)


for line in body: print(line)
print ('</svg>')

"""


