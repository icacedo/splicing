import isomod as im
import sys

fasta = sys.argv[1]
gff = sys.argv[2]
exons = sys.argv[3]

seqid, seq = im.read_fasta(fasta)

print(seq)

dons, accs = im.read_gff_sites(seq, gff, gtag=True)

print(dons, accs)

re_seq = im.read_txt_seqs(exons)

#print(re_seq)

exinlens, a, b, g, size_limit = im.fdist_params(re_seq)

scores, y_values = im.memoize_fdist(data, a, b, g, 500)

print(len(values))
print(max(values))
print(min(values))
print(max(exinlens), min(exinlens))







