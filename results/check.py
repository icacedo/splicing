import isoform
seq = 'ATATATCGTCGATCAGCCGTATATATCGTCGATCAGCCGT'
print(len(seq))
accs, dons = isoform.gtag_sites(seq, 0, 0)
print(accs, dons)
