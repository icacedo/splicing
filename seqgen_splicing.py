import random

random.seed(1)

# generate n random sequences of length l
# specify number of GT and AG sites to be randomly included
# 0.5% chance that a position starts with GT or AG
# AT/GC content is 40/60

# what to do if GT and AG sites are randomly generated?

def ranseq_splicing(seq_l,seq_n,don_n,acc_n):
	for i in range(seq_l):
		r=round(random.random(),2)
		don_or_acc=random.randint(1,2)
		if r<=0.20 and don_or_acc==1: 
			print('GT')
		if r<=0.20 and don_or_acc==2:
			print('AG')
		
	
seq_n=10
seq_l=100	
don_n=30
acc_n=20
ranseq_splicing(seq_l,seq_n,don_n,acc_n)
