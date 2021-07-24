import argparse
import itertools as it
import seqlib as sl

parser=argparse.ArgumentParser(description='splice this sequence')

parser.add_argument('--fasta', required=True, type=str, 
	help='fasta file to splice')

arg=parser.parse_args()

#print(arg.fasta)

don_acc_pos={}
gt=0
ag=0
for i,j in sl.read_fasta(arg.fasta):
	for p in range(len(j)):
		if j[p:p+2]=='GT':
			don_acc_pos[p]='GT'
			gt+=1
		if j[p:p+2]=='AG':
			don_acc_pos[p]='AG'
			ag+=1
		else:continue
#print(don_acc_pos, len(don_acc_pos), gt+ag)
'''
for i,j in don_acc_pos.items():
	print(i,j)
	print(don_acc_pos[i])
	print(i)
	print(don_acc_pos)
	break
	
pos = don_acc_pos.keys()
print(pos)
type = don_acc_pos.values()
print(type)
'''
# start with this loop after to include all lengths of tuples
# replace 2 with i in second loop
#for i in range(len(don_acc_pos)):
for j in it.combinations(don_acc_pos.keys(), 2):
	#print(j)
	for k in j:
		#print(don_acc_pos[k])
		if don_acc_pos[k]=='AG':print('############')
		else: print('@@@@@@@@@')
		break

	
		
