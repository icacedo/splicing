import isoform2

seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'

maxs = 3
minin = 3
minex = 3
flank = 3

dons, accs = isoform2.gtag_sites(seq, flank, minex)

print(dons, accs)

introns = []
for d in dons:
    for a in accs:
        introns.append((d, a))

print(introns)