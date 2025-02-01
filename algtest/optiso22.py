import argparse
import subprocess
import random
import sys
import json

# ln -s isoform.py into working directory
import isoform2
from isoform2 import Locus
'''
parser = argparse.ArgumentParser()
parser.add_argument('config', type=str, metavar='<json>',
    help='configuration file with all genes')
#parser.add_argument('apc_gen_gene', type=str, metavar='<file>',
#    help='apc generated gff file')
parser.add_argument('--apc_prog', required=False, type=str, default='geniso2',
	metavar='<exec>', help='path to apc program [%(default)s]')
parser.add_argument('--cmp_prog', required=False, type=str, default='cmpiso',
    metavar='<exec>', help='path to isoform comparison program \
        [%(default)s]')

args = parser.parse_args()


with open('config_optiso2.json', 'r') as file:
    data = json.load(file)
'''

seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'

maxs = 3
minin = 3
minex = 3
flank = 3

dons, accs = isoform2.gtag_sites(seq, flank, minex)

print('gtag sites:', dons, accs)

# try all permutations of acceptors, given a single donor site
# so given donor 6, 
# [6, 17, 22] [13, 28]
# all possible introns
# [6, 13] 1
# [6, 28] 2
# [17, 13] 3 invalid
# [17, 28] 4
# [22, 13] 5 invalid
# [22, 28] 6
# all possible combinations of introns (each has unique index)
# 1, 2 invalid
# 1, 3
# 1, 4
# 1, 5
# and so on...


import sys

#fasta = sys.argv[1]

introns = []
for d in dons:
    for a in accs:
        if d >= a: continue
        introns.append((d, a))

print('introns:', introns)

# donor acceptor pairs
# [6, 17, 22] [13, 28]
# all possible introns

# [6, 13]
# [6, 28]
# [6, 13] [17, 13] !
# [6, 13] [17, 28]
# [6, 13] [17, 28] [22, 13] !
# [6, 13] [17, 28] [22, 28] !
# [17, 13] !
# [17, 28]
# [22, 28]
# [6, 13] [22, 28]

import copy

print('##########')



# from chatgpt
# but it's incorrect
# missing a few possible isoforms
'''
def build_isoforms(dons, accs, introns, imin, emin, limit=None, countonly=False, isocount=0):
    if countonly and limit and isocount >= limit:
        return
    
    don = dons[0]
    for aix, acc in enumerate(accs):
        if acc - don + 1 < imin:
            continue
        intron = (don, acc)

        # legit isoform as is, print it
        if countonly and limit and isocount >= limit:
            return
        iso = copy.copy(introns)
        iso.append(intron)
        print(iso)  # Print isoform

        # also extend it
        descendable = False
        for dix, ndon in enumerate(dons):
            elen = ndon - acc - 1
            if elen >= emin:
                ext = copy.copy(iso)
                build_isoforms(dons[dix:], accs[aix:], ext, imin, emin, limit, countonly, isocount)

build_isoforms(dons, accs, [], 1, 1)
'''

# works on test seq
'''
model = sys.argv[1]
model = isoform2.read_splicemodel(model)

weights = {
    'wacc': 1.0,
    'wdon': 1.0,
    'wexs': 1.0,
    'wins': 1.0,
    'wexl': 1.0,
    'winl': 1.0,
    'winf': 1.0
}

constraints = {
    'min_intron': 1,
    'min_exon': 1,
    'flank': 1
}


locus = Locus('test', seq, model, constraints, weights)

locus.gff(sys.stdout)
'''

# [6, 17, 22] [13, 28]
'''
[6, 13]
[6, 13] [17, 28]
[6, 28]
'''

dons = [6, 17, 22]
accs = [13, 28]

for aix, acc in enumerate(accs):
    intron = (dons[0], acc)
    print(intron)
    for dix, don in enumerate(dons):
        intron = (don, acc)

'''
[6, 13]
[6, 28]
[6, 13] [17, 13] 
[6, 13] [17, 28]
[6, 13] [17, 28] [22, 13] 
[6, 13] [17, 28] [22, 28] 
[17, 13] 
[17, 28]
[22, 28]
[6, 13] [22, 28]
'''