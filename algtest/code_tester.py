# testing some of ian's code so i can write about

# calculating the inttron cost
# isoforms/modelbuilder

import glob
import argparse
from grimoire.genome import Reader
import math

parser = argparse.ArgumentParser()
parser.add_argument('dir')

args = parser.parse_args()

exons = []
introns = []
exonsum = 0
for ff in glob.glob(f'{args.dir}/*.fa'):
    gf = ff[:-2] + 'gff3'
    genome = Reader(gff=gf, fasta=ff)
    tx = next(genome).ftable.build_genes()[0].transcripts()[0]
    for f in tx.exons:
        exons.append(f.seq_str())
        exonsum += f.end - f.beg + 1
    for f in tx.introns:
        iseq = f.seq_str()
        introns.append(iseq)

# same inf as in worm.splicemodel
print(len(introns), exonsum)
inf = len(introns)/exonsum
print(math.log2(inf/0.25))
