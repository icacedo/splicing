import random
import seqlib as sl
import numpy as np

# n is the length of an individual sequence
# N is the number of sequences
def simple_seqgen(n, N):

	seq_list = []
	for i in range(1, N + 1):
		#id='>s'+str(i)
		seq=sl.random_dna(n,0.20,0.30,0.30,0.20)
		seq_list.append(seq)
	return(seq_list)

def make_pfm(seqs, n, N):

	pfm = [
		{'A':0, 'C':0, 'G':0, 'T':0} for i in range(n)
		]
	for seq in seqs:
		position = 0
		for nt in seq:
			pfm[position][nt] += 1
			position += 1
	for i in range(len(pfm)):
		for nt in pfm[i]:
			pfm[i][nt] = pfm[i][nt]/N
	return(pfm)

def make_pwm(seqs, n, N, pseudo_count):

	pfm_arr = np.zeros(shape=(4, n))
	A = 0
	C = 0
	G = 0
	T = 0
	for seq in seqs:
		position = 0
		for nt in seq:
			if nt.upper() == 'A':
				pfm_arr[0][position] += 1
				A += 1
			if nt.upper() == 'C':
				pfm_arr[1][position] += 1
				C += 1
			if nt.upper() == 'G':
				pfm_arr[2][position] += 1
				G += 1
			if nt.upper() == 'T':
				pfm_arr[3][position] += 1
				T += 1
			position += 1
	pfm_arr = np.divide(pfm_arr, N) + pseudo_count

	background = A + C + G + T
	Af = A / background
	Cf = C / background
	Gf = G / background
	Tf = T / background

	pwm_arr = np.zeros(shape=(4, n))
	row = 0
	for x in pfm_arr:
		col = 0
		for y in x:
			if row == 0: 
				pwm_arr[row][col] = np.log2(pfm_arr[row][col] / Af)
			if row == 1:
				pwm_arr[row][col] = np.log2(pfm_arr[row][col] / Cf)
			if row == 2:
				pwm_arr[row][col] = np.log2(pfm_arr[row][col] / Gf)
			if row == 3:
				pwm_arr[row][col] = np.log2(pfm_arr[row][col] / Tf)
			col += 1
		row += 1
	return(pwm_arr)
	






















	
