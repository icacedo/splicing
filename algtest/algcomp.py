import argparse
import isoform2
import copy

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
                    help='single APC gene to test')
parser.add_argument('--model')

args = parser.parse_args()

name, seq = next(isoform2.read_fasta(args.fasta))

flank = 99
emin = 25
imin = 35

'''
seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'

maxs = 3
imin = 3
emin = 3
flank = 3
'''
dons, accs = isoform2.gtag_sites(seq, flank, emin)

# Ian's algorithm

isos = []
introns = []
def buildIsos(dons, accs, introns):

    don = dons[0]
    for aix, acc in enumerate(accs):
        if acc - don + 1 < imin: continue
        intron = (don, acc)
        iso = copy.copy(introns)
        iso.append(intron)
        print(iso, '***')
        isos.append(iso)

        for dix, ndon in enumerate(dons):
            elen = ndon - acc - 1
            if elen >= emin:
                ext = copy.copy(iso)
                buildIsos(dons[dix:], accs[aix:], ext)

buildIsos(dons, accs, introns)

print(isos, '*')

print('next section')
# recursive backtracking/leetcode #78

introns = []
for d in dons:
    for a in accs:
        ilen = a-d+1
        if ilen >= imin:
            introns.append((d, a))

def solCheck(sol):

    if sol == []: 
        return 
    
    sites = []
    for tup in sol:
        sites.append(tup[0])
        sites.append(tup[1])

    if len(sites) != len(set(sites)):
        return 
    
    else:
        return sol

def subsets(introns):
    n = len(introns)
    res, sol = [], []

    def backtrack(i):
        if i == n:
            res.append(sol[:])
            return

        backtrack(i+1)

        sol.append(introns[i])
        backtrack(i+1)
        sol.pop()

    backtrack(0)
    return res

res = subsets(introns)

'''
print('#####')

# works...but how to code in imin/emin/flank?
for r in res:
    if solCheck(r): 
        print(r)

print('#####')

name = 'test'
seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'

maxs = 3
minin = 3
minex = 3
flank = 3

dons, accs = isoform2.gtag_sites(seq, flank, minex)


from isoform2 import Locus

model = isoform2.read_splicemodel(args.model)

constraints = {
    'min_intron': 3,
    'min_exon': 3,
    'flank': 3
}

weights = {
	'wacc': 1.0,
	'wdon': 1.0,
	'wexs': 1.0,
	'wins': 1.0,
	'wexl': 1.0,
	'winl': 1.0,
	'winf': 1.0,
}

locus = Locus(name, seq, model, constraints, weights)

import sys
#locus.gff(sys.stdout)


(18, 29)
(23, 29)
(7, 29)
(7, 14)
(7, 14) (18, 29)
(7, 14) (23, 29)


# backtracking is way slower than ian's method

 

'''

