import seqlib 
import random
<<<<<<< HEAD
import numpy
=======
import isoform
import numpy as np
>>>>>>> 1c52f8e4e843dd6ff0f258f09309ab64284fa1fa

random.seed(2)

def simple_seqgen(n,N):
	seq_list=[]
	for i in range(1,N+1):
		#id='>s'+str(i)
		seq=seqlib.random_dna(n,0.20,0.30,0.30,0.20)
		seq_list.append(seq)
	return(seq_list)

n=10
N=5

seqs=simple_seqgen(10,5)

def make_pwm(seqs,n):
	pwm=[{} for i in range(n)]
	for seq in seqs:
		position=0
		for nt in seq:
			if nt not in pwm[position]:
				pwm[position][nt]=1
			else:
				pwm[position][nt]+=1
			position+=1
	for i in range(len(pwm)):
		for nt in pwm[i]:
			pwm[i][nt]=pwm[i][nt]/N
	return(pwm)

#print(make_pwm(seqs,n))

pfm_arr=np.zeros(shape=(4,n))
A=0
C=0
G=0
T=0
for seq in seqs:
	position=0
	for nt in seq:
		if nt.upper()=='A':
			pfm_arr[0][position]+=1
			A+=1
		if nt.upper()=='C':
			pfm_arr[1][position]+=1
			C+=1
		if nt.upper()=='G':
			pfm_arr[2][position]+=1
			G+=1
		if nt.upper()=='T':
			pfm_arr[3][position]+=1
			T+=1
		#else: 
		#	print('unknown character')
		#	sys.exit()
		position+=1
pfm_arr=np.divide(pfm_arr,N)
print(pfm_arr)

background=A+C+G+T
Af=A/background
Cf=C/background
Gf=G/background
Tf=T/background

pwm_arr=np.zeros(shape=(4,n))
pseudo=0.001
r=0
for x in pfm_arr:
	c=0
	for y in x:
		if r==0:
			if pfm_arr[r][c]==0:
				pfm_arr[r][c]+=0.001
			pwm_arr[r][c]=np.log2(pfm_arr[r][c]/Af)
		if r==1:
			if pfm_arr[r][c]==0:
				pfm_arr[r][c]+=0.001
			pwm_arr[r][c]=np.log2(pfm_arr[r][c]/Cf)
		if r==2:
			if pfm_arr[r][c]==0:
				pfm_arr[r][c]+=0.001
			pwm_arr[r][c]=np.log2(pfm_arr[r][c]/Gf)
		if r==3:
			if pfm_arr[r][c]==0:
				pfm_arr[r][c]+=0.001
			pwm_arr[r][c]=np.log2(pfm_arr[r][c]/Tf)
		c+=1
	r+=1
print(pwm_arr)
























	
