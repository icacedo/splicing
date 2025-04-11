import argparse
import isoform2
import copy

parser = argparse.ArgumentParser()
parser.add_argument('fasta')

args = parser.parse_args()

name, seq = next(isoform2.read_fasta(args.fasta))

flank = 99
emin = 25
imin = 35

seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'

maxs = 3
imin = 3
emin = 3
flank = 3

dons, accs = isoform2.gtag_sites(seq, flank, emin)

isos = []
def buildIsos(dons, accs, introns):

	don = dons[0]
	for aix, acc in enumerate(accs):
		if acc - don + 1 < imin: continue
		intron = (don, acc)
		iso = introns + [intron]
		isos.append(iso)
		for dix, ndon in enumerate(dons):
			elen = ndon - acc - 1
			if elen >= emin:
				ext = iso
				buildIsos(dons[dix:], accs[aix:], ext)
print(dons)
print(accs)
buildIsos(dons, accs, [])
print(isos)


# practice recursion here
# https://www.geeksforgeeks.org/backtracking-algorithms/

nums = [1, 2, 3]

subset = []
res = []
def subsetRecur(i, arr, res, subset):

	if i == len(arr):
		res.append(list(subset))
		return

	subset.append(arr[i])
	subsetRecur(i + 1, arr, res, subset)

	subset.pop()
	subsetRecur(i + 1, arr, res, subset)

subsetRecur(0, nums, res, subset)

print(res)





