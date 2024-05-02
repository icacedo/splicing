import isomod as im
import sys
import apc_model_lib as aml

fasta = sys.argv[1]
gff = sys.argv[2]
wb_dir = sys.argv[3]

seqid, seq = im.read_fasta(fasta)

dons, accs = im.read_gff_sites(seq, gff)

#all_es, all_is, all_ds, all_as = im.get_all_tn_seqs(wb_dir)

#im.get_top_exins(seq, gff)


'''
exlens, a, b, g, size_limit = im.fdist_params(all_es)

scores, y_values = im.memoize_fdist(exlens, a, b, g, 1000)
print(exlens)
print(len(y_values))
print(max(y_values))
print(min(y_values))
print(max(exlens), min(exlens))

new_ys = im.zero_prob(exlens, y_values)

#print(new_ys)
totl = 0
for i in y_values:
    totl += i
print(totl)

tot = 0
for i in new_ys:
    tot += i
print(tot)

elens = []
for seq in all_es:
    if len(seq) <= 20:
        print(seq)
'''
# isoforms
# 3a.1
# 3c.1
# 3c.2









