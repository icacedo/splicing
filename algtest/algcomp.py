import argparse
import isoform2
import copy

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
                    help='single APC gene to test')

args = parser.parse_args()
'''
name, seq = isoform2.read_fasta(args.fasta)

flank = 99
minex = 25
minin = 35
'''
seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'

maxs = 3
imin = 3
emin = 3
flank = 3

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
        isos.append(iso)

        for dix, ndon in enumerate(dons):
            elen = ndon - acc - 1
            if elen >= emin:
                ext = copy.copy(iso)
                buildIsos(dons[dix:], accs[aix:], ext)

buildIsos(dons, accs, introns)

print(isos)


print('#####')
# recursive backtracking/leetcode #78

introns = []
for d in dons:
    for a in accs:
        ilen = a-d+1
        if ilen >= imin:
            introns.append((d, a))

def solCheck(sol, emin, imin, flank):

    if sol == []: 
        return False

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

print(subsets(introns))

print(introns)
print('#####')

valid = []
flank = 3
for icoors in introns:
    if icoors[0] > flank and icoors[1] <= len(seq) - flank + 1:
        valid.append(icoors)
print(valid)

def dunno(valid):
    seenD = []
    seenA = []
    possible = []
    for intron in valid:
        if intron[0] not in seenD:
            seenD.append(intron[0])
            if intron[1] not in seenA:
                seenA.append(intron[1])
                possible.append(intron)
    print(possible)

dunno(valid)

dunno(valid[1:])

subtest1 = [(6, 13), (6, 28), (17, 28)]
subtest2 = [(6, 13), (17, 28)]

def uniqueSites(sublist):

    sites = []
    for item in sublist:
        sites.append(item[0])
        sites.append(item[1])

    if len(sites) != len(set(sites)): return False
    else: return True

print(uniqueSites(subtest2))
print(uniqueSites(subtest1))

print('#################')
def sub(a):
    sublists = []
    for i in range(len(a)):
        for j in range(i + 1, len(a) + 1):  
            print(a[i:j], '$$$')
            print(uniqueSites(a[i:j]))
            sublists.append(a[i:j])
    return sublists

print(sub(introns))

print('############')

print(introns)

sublists = []
for i in range(len(introns)):
    for j in range(i + 1, len(introns) + 1):
        sublists.append(introns[i:j])

for s in sublists:
    print(s)
    if uniqueSites(s):
        print(s)





