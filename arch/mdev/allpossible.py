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

# convert output of combinations() to a dictionary
def convert(tup_don,tup_acc):
	combos={}
	for p in tup_don:
		combos[p]='GT'
	for p in tup_acc:
		combos[p]='AG'
	yield combos

'''
for i in range(len(don)):
	for j in it.combinations(don, 2):
		for k in it.combinations(acc, 2):
			for c in convert(j,k):
				c_list=sorted(c)
				key_first=c_list[0]
				key_last=c_list[2+2]
				# filter out combos starting with an acceptor site
				if c[key_first]=='AG': continue
				if c[key_last]=='GT': continue
				break
'''
total=0
for d in range(len(don)):
	for dsites in it.combinations(don, d):
		for asites in it.combinations(acc, d):
			for c in convert(dsites,asites):
				sorted_positions=sorted(c)	
				splice_sites={}
				prev=''
				for p in range(len(sorted_positions)):
					if c[sorted_positions[p]]=='AG' and prev=='': 
						prev='AG'
					if c[sorted_positions[p]]=='GT' and prev=='': 
						splice_sites[sorted_positions[p]]='GT'
						prev='GT'
					if c[sorted_positions[p]]=='GT' and prev=='AG': 
						splice_sites[sorted_positions[p]]='GT'
						prev='GT'
					if c[sorted_positions[p]]=='AG' and prev=='GT': 
						splice_sites[sorted_positions[p]]='AG'
						prev='AG'
					if c[sorted_positions[p]]=='GT' and prev=='GT':
						prev='GT'
					if c[sorted_positions[p]]=='AG' and prev=='AG':
						prev='AG'
					if len(splice_sites)>=2:
						# there are a lot of duplicates
						#print(splice_sites)
						total+=1
print(total)

			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
				
			

	
		
