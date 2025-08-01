# practice coding the Smith-Waterman algorithm
# quadratic memory-allocated the matrix into memory
# as opposed to keeping only the two rows you are working on (Meyers-Miller)
# modify with ideas from posted paper 
# look for ian's blast primer, not in the lab repo

seq1 = 'TGTTACGG'
seq2 = 'GGTTGACTA'

seq1 = 'CGTGAATTCAT'
seq2 = 'GACTTAC'

matrix = [[0 for x in range(len(seq1)+1)] for x in range(len(seq2)+1)]
traces = [[None for x in range(len(seq1)+1)] for x in range(len(seq2)+1)]

match = 3
mismatch = -3 
gap = -2

# get highest score
h = [0, 0, 0]
for i in range(1, len(matrix)):
	for j in range(1, len(matrix[i])):

		diag = matrix[i-1][j-1]
		left = matrix[i][j-1]
		up = matrix[i-1][j] 
		
		diag_score = 0
		if seq1[j-1] == seq2[i-1]: 
			diag_score += match + diag
		else:
			diag_score += mismatch + diag
			
		if diag_score > h[0]:
				h = [diag_score, i, j]
			
		left_score = left + gap
		up_score = up + gap
		
		if diag_score >= left_score and diag_score >= up_score and diag_score > 0:
			matrix[i][j] = diag_score
			traces[i][j] = 'diag'
		elif left_score >= up_score and left_score > 0:
			matrix[i][j] = left_score
			traces[i][j] = 'left'
		elif up_score > 0:
			matrix[i][j] = up_score
			traces[i][j] = 'up'
		else:
			matrix[i][j] = 0

align1 = []
align2 = []
while matrix[h[1]][h[2]] > 0:
	#print(h, traces[h[1]][h[2]], seq1[h[2]-1], seq2[h[1]-1])
	if traces[h[1]][h[2]] == 'diag':
		align1.append(seq1[h[2]-1])
		align2.append(seq2[h[1]-1])
		ni = h[1]-1
		nj = h[2]-1
		h = [matrix[ni][nj], ni, nj]
		continue
	if traces[h[1]][h[2]] == 'up':
		align1.append('-')
		align2.append(seq2[h[1]-1])
		ni = h[1]-1
		nj = h[2]
		h = [matrix[ni][nj], ni, nj]
		continue
	if traces[h[1]][h[2]] == 'left':
		align1.append(seq1[h[2]-1])
		align2.append('-')
		ni = h[1]
		nj = h[2]-1
		h = [matrix[ni][nj], ni, nj]
		continue
		
align1.reverse()
align2.reverse()
print(''.join(align1))
print(''.join(align2))
		
# print matrix 
def print_matrix(matrix):
	
	score_lengths = []
	for row in matrix:
		for score in row:
			score_lengths.append(len(str(score)))

	width = max(score_lengths)

	# dynamically change width of each cell displayed
	def add_space(string, width):
		
		string_len = len(string)
		while string_len < width:
			string += ' '
			string_len += 1
			
		return string

	seq_string = '  ' + (' ' * (width+1))
	for nuc in seq1:
		seq_string += add_space(f'{nuc}', width+1)

	print(seq_string)
	for i, row in enumerate(matrix):
		if i == 0:
			row_string = '  '
		else:
			row_string = f'{seq2[i-1]} '
		for score in row:
			row_string += add_space(f'{score}', width+1)
		print(row_string)


print_matrix(matrix)

print_matrix(traces)

	

	
	
	
	
	
	
