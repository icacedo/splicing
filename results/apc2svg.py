#!/usr/bin/env python3

import argparse
import subprocess
import sys

apc = 'build/apc'

def isoforms(dna, prog, fasta, mdir, don, acc, emm, imm, elen, ilen, icost):
	stuff = subprocess.run([prog, fasta,
		'--dpwm', f'{mdir}/donor.pwm',    '--wdpwm', don,
		'--apwm', f'{mdir}/acceptor.pwm', '--wapwm', acc,
		'--emm',  f'{mdir}/exon.mm',      '--wemm',  emm,
		'--imm',  f'{mdir}/intron.mm',    '--wimm',  imm,
		'--elen', f'{mdir}/exon.len',     '--welen', elen,
		'--ilen', f'{mdir}/intron.len',   '--wilen', ilen,
		'--icost', str(icost)], capture_output=True, text=True).stdout

with open('results/arg.tsv') as fp:
	header = fp.readline()
	for line in fp.readlines():
		seq, don, acc, emm, imm, elen, ilen, icost, fit = line.split()
		genome = grimoire.genome.Reader(
			fasta=f'{apc}/{seq}.fa',
			gff=f'{apc}/{seq}.gff3')
		chrom = next(genome)
		gene = chrom.ftable.build_genes()[0]
		tx = gene.transcripts()[0]
		

completely unfinished

