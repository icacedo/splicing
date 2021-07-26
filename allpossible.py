import argparse
import itertools as it
import seqlib as sl

parser=argparse.ArgumentParser(description='splice this sequence')

parser.add_argument('--fasta', required=True, type=str, 
	help='fasta file to splice')

arg=parser.parse_args()

don=[]
acc=[]
for i,j in sl.read_fasta(arg.fasta):
	for p in range(len(j)):
		if j[p:p+2]=='GT':
			don.append(p)
		if j[p:p+2]=='AG':
			acc.append(p)
		else:continue

def convert(tup_don,tup_acc):
	#don={}
	#acc={}
	combos={}
	for p in tup_don:
		#don[p]='GT'
		combos[p]='GT'
	for p in tup_acc:
		#acc[p]='AG'
		combos[p]='AC'
	#yield don,acc
	yield combos
	#### need to filter out libraries with UNPAIRED GT/AG sites
for i in range(len(don)):
	for j in it.combinations(don, 2):
		for k in it.combinations(acc, 2):
			for c in convert(j,k):
				c_list=sorted(c)
				key_first=c_list[0]
				key_last=c_list[2+2]
				# filter out combos starting with an acceptor site
				if c[key_first]=='AG': continue
				if c[key_last]=='GC': continue
				break
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
				
			

	
		
