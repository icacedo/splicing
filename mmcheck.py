import sys
import modelib as ml
import isoform_fixed as isof

seq = 'AACATGACCGTTGCGAGCTACCGTCACATTAGCTCGGAGCCCTATATA'

# intron mm file
fp1 = sys.argv[1]

# donor pwm file
fp2 = sys.argv[2]

# acceptor pwm file
fp3 = sys.argv[3]

dons, accs = ml.get_gtag(seq)

mm_probs, mm_scores = ml.read_exin_mm(fp1)

dppm, dpwm = ml.read_pwm(fp2)

appm, apwm = ml.read_pwm(fp3)

apc_isoforms, trials = ml.apc(dons, accs, 100, 3, 4, 5, seq)

for iso in apc_isoforms:
	exon_seqs, intron_seqs = ml.get_exin_seqs(iso, seq)
	break

print(len(dpwm), len(apwm))

len_dpwm = 2
len_apwm = 2


intron_score_total = 0
for inseq in intron_seqs:
	intron_score = 0
	inseq = inseq[len_dpwm:-len_apwm]
	for i in range(len(inseq)):
		if len(inseq[i:i+k]) == k:
			kmer = 
	




#def get_mm_score(exin_seqs, exin_mm, dpwm=None, apwm=None):
	
	#if dpwm and apwm:
		
		
