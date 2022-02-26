import modelib
import seqlib
import sys

# coding a markov chain model
# starting with viterbi? 

'''
# code to make simple test file
fp = open('test_seqs.txt', 'w')
seqs = modelib.simple_seqgen(10,5)

for i in range(len(seqs)):
	fp.write(seqs[i])
	fp.write('\n')
'''
'''
# maybe don't need this for viterbi
# test files are in the viterbi repo
fp = open(sys.argv[1])

for line in fp.readlines():
	print(line.rstrip())
'''

##### viterbi training section ################################################
exons = sys.argv[1]
introns = sys.argv[2]

# this will get global nucleotide frequencies
# make_pfm in modelib.py gets positional frequencies
def nuc_freqs(fasta):
	nfreqs = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
	total_lens = 0
	seq_count = 0

	for label, seq in seqlib.read_fasta(fasta):
		total_lens += len(seq)
		seq_count += 1
		for n in seq:
			nfreqs[n] += 1
		

	total_n = 0
	for i in nfreqs:
		total_n += nfreqs[i]
	for i in nfreqs:
		nfreqs[i] = nfreqs[i]/total_n
	
	avg_len = total_lens/seq_count
	
	return nfreqs, avg_len

print(nuc_freqs(exons))
print(nuc_freqs(introns))


