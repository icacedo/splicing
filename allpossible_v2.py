import argparse
import itertools as it
import seqlib as sl

parser=argparse.ArgumentParser(description='splice this sequence')
parser.add_argument('--fasta', required=False, type=str, 
	help='fasta file to splice')
parser.add_argument('--minexon', required=False, type=int, 
	default=25, help='minimum exon length')
parser.add_argument('--minintron', required=False, type=int,
	default=35, help='minimum intron length')
arg=parser.parse_args()

fasta=arg.fasta
minexon=arg.minexon
minintron=arg.minintron

def find_sites_fasta(seq,minexon):

	don_sites=[]
	acc_sites=[]
	for i,j in sl.read_fasta(seq):
		# range needs to take into account minimum exon size
		for p in range(minexon,len(j)-minexon):
			if j[p:p+2]=='GT':
				don_sites.append(p)
			if j[p:p+2]=='AG':
				# +1 added to get to the end of the acceptor sequence
				acc_sites.append(p+1)
			else:continue
		yield i,don_sites,acc_sites
		don_sites=[]
		acc_sites=[]
	
def find_sites(seq,minexon):
	
	don_sites=[]
	acc_sites=[]
	for i in range(minexon,len(seq)-minexon):
		if seq[i:i+2]=='GT':
			don_sites.append(i)
		if seq[i:i+2]=='AG':
			acc_sites.append(i)
		else:continue
	yield don_sites,acc_sites

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

def all_possible(don_sites,acc_sites,minexon,minintron):

	info={
	'total_isoforms':0,
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
				info['total_isoforms']+=1
	return info

def ranseq_test(min_seq,max_seq,seq_step,n_seqs):
		
	for length in range(min_seq,max_seq,seq_step):
		for amount in range(n_seqs):
			# 60% GC content
			yield sl.random_dna(length,0.20,0.30,0.30,0.20)
# use this chunk to test using many random sequences
'''
id=0
for seq in ranseq_test(100,300,50,10):
	for d,a in find_sites(seq,minexon):
		id+=1
		print(id,len(seq),all_possible(d,a,minexon,minintron))
'''

# use this chunk to test using a fasta file
for i,d,a in find_sites_fasta(fasta,minexon):
	print(i,all_possible(d,a,minexon,minintron))


					
				
				
				

