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
				ext = iso
				buildIsos(dons[dix:], accs[aix:], ext)
print(dons)
print(accs)
buildIsos(dons, accs, introns)
print(isos)

