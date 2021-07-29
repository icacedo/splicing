# rewrite/methods check of code using combinations()
# need to integrate the following parameters:
# seq, minin, minex, dons, accs
# output should be the same as isoform.py/isoformer.py

#import argparse
import itertools as it
import seqlib as sl

'''
parser=argparse.ArgumentParser(description='splice this sequence')
parser.add_argument('--fasta', required=True, type=str, 
	help='fasta file to splice')
arg=parser.parse_args()

# all sites need to be sorted by order of position
# this method naturally sorts sites
don_sites=[]
acc_sites=[]
for i,j in sl.read_fasta(arg.fasta):
	for p in range(len(j)):
		if j[p:p+2]=='GT':
			don_sites.append(p)
		if j[p:p+2]=='AG':
			acc_sites.append(p)
		else:continue
'''

def makesnosense(dons,accs):

	for d,a in zip(dons,accs):
		if d>a: return True
	return False
	
def short_introns(dons,accs,minintron):

	for d,a in zip(dons,accs):
		intron_length=a-d
		if intron_length<minintron:
			return True
	return False
	
def short_exons(dons,accs,minexon):

	for i in range(1,len(dons)):
		exon_begin=accs[i-1]+1
		exon_end=dons[i]-1
		exon_length=exon_end-exon_begin
		if exon_length<minexon: 
			return True
	return False

def all_possible(don_sites,acc_sites,minintron,minexon):

	info={
	'trials':0,
	'n_dsites':len(don_sites),
	'n_asites':len(acc_sites),
	'coor_fails':0,
	'intron_fails':0,
	'exon_fails':0
	}	

	for n in range(1,len(don_sites)+1):
		for d in it.combinations(don_sites,n):
			for a in it.combinations(acc_sites,n):
				info['trials']+=1
				if makesnosense(d,a): 
					info['coor_fails']+=1
					continue	
				if short_introns(d,a,minintron):
					info['intron_fails']+=1
					continue	
				if short_exons(d,a,minexon):
					info['exon_fails']+=1
					continue				
				#print(d,a)
	return info

def ranseq_test(min_seq,max_seq,seq_step,n_seqs):
		
	for length in range(min_seq,max_seq,seq_step):
		for amount in range(n_seqs):
			# 60% GC content
			yield sl.random_dna(length,0.20,0.30,0.30,0.20)

def find_sites(seq):
	
	don_sites=[]
	acc_sites=[]
	for i in range(len(seq)):
		if seq[i:i+2]=='GT':
			don_sites.append(i)
		if seq[i:i+2]=='AG':
			acc_sites.append(i)
		else:continue
	yield don_sites,acc_sites
	
for seq in ranseq_test(100,300,50,10):
	for d,a in find_sites(seq):
		print(all_possible(d,a,20,40))


					
				
				
				

