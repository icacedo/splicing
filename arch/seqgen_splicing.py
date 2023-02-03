import random

random.seed(1)

# generate n random sequences of length l
# specify number of GT and AG sites to be randomly included
# 0.5% chance that a position starts with GT or AG
# AT/GC content is ignored, too hard for now

# what to do if GT and AG sites are randomly generated?

def ranseq_splicing(seq_l,seq_n,don_n,acc_n):
	for i in range(seq_l):
		r=round(random.random(),2)
		print(r)
		don_or_acc=random.randint(1,2)
		if r<=0.20 and don_or_acc==1: 
			print('GT')
		if r<=0.20 and don_or_acc==2:
			print('AG')
		if r>0.20:
			if r>=0.21 and r<=0.40:
				print('A')
			if r>=0.41 and r<=0.60:
				print('C')
			if r>=0.61 and r<=0.80:
				print('G')
			# why does the output change if I don't specify r>=0.81?
			if r>=0.81:
				print('T')
			
		
	
seq_n=10
seq_l=100	
don_n=30
acc_n=20
ranseq_splicing(seq_l,seq_n,don_n,acc_n)
