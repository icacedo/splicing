import argparse
import itertools as it
import seqlib as sl

parser=argparse.ArgumentParser(description='splice this sequence')

parser.add_argument('--fasta', required=True, type=str, 
	help='fasta file to splice')

arg=parser.parse_args()

print(arg.fasta)

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
print(len(don_acc_pos), gt+ag)
		
