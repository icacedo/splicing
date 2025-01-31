import argparse
import subprocess
import random
import sys
import json

# ln -s isoform.py into working directory
import isoform
from isoform import Locus
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

dn, ac = isoform.gtag_sites(seq, flank, minex)

print(dn, ac)

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
for d in dn:
    for a in ac:
        if d >= a: continue
        introns.append((d, a))



def perm(list, start, end):
    if (start == end):
        print(list)
    else:
        for i in range(start, end + 1):
            list[start], list[i] = list[i], list[start]
            perm(list, start + 1, end)
            list[start], list[i] = list[i], list[start]

perm(introns, 0, 2)


