import modelib as ml
import sys
import isoform_fixed as isoform
import gzip

fasta = sys.argv[1]
len_model = sys.argv[2]
exon_seqs = sys.argv[3]

def get_seqs(filename):
	seqs = []
	with gzip.open(filename, 'rt') as fp:
		for seq in fp.readlines():
			seq = seq.rstrip()
			seqs.append(seq.rstrip())
	return seqs

exons = get_seqs(exon_seqs)

elen = isoform.create_len(exons, 15, 500)

model_read = isoform.read_len(len_model)
print(model_read)

# value from isoform.read_len
tail = isoform.find_tail(0.00035, 500)
print(tail)

print(len(elen))

for seqid, seq in ml.read_fastas(fasta):
	seqid = seqid
	seq = seq

dons, accs = ml.get_gtag(seq)


'''
apc_isoforms, trials = ml.apc(dons, accs, 3, 25, 25, 100, seq)

re_elen_pdf, re_elen_log2 = ml.read_exin_len(len_model)

for iso in apc_isoforms:
	exon_lengths, intron_lengths = ml.get_exin_lengths(iso)
	elen_score = ml.get_len_score(exon_lengths, re_elen_log2)
	print(elen_score)
	print(iso)
'''







