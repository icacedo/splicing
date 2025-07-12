# practice coding the Smith-Waterman algorithm
# quadratic memory-allocated the matrix into memory
# as opposed to keeping only the two rows you are working on (Meyers-Miller)
# modify with ideas from posted paper 
# look for ian's blast primer, not in the lab repo

seq1 = 'GAATTCAGTTA'
seq2 = 'GGATCGA'

matrix = [[0 for x in range(len(seq1)+1)] for x in range(len(seq2)+1)]

match = 5
mismatch = -3 
gap = -4

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
			
		left_score = left + gap
		up_score = up + gap
		
		score = max(diag_score, left_score, up_score)
		if score < 0:
			matrix[i][j] = 0
		else:
			matrix[i][j] = score

# print matrix 
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

# get tracebacks
highest = {}
for i in range(len(matrix)):
	print(max(matrix[i]), i)
	highest[i] = max(matrix[i])
	

for item in highest.items():
	print(item)
	
	
	
	
	
	
	
	
