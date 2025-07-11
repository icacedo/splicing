# practice coding the Smith-Waterman algorithm
# quadratic memory-allocated the matrix into memory
# as opposed to keeping only the two rows you are working on (Meyers-Miller)
# modify with ideas from posted paper 
# look for ian's blast primer, not in the lab repo

seq1 = 'GAATTCAGTTA'
seq2 = 'GGATCGA'

seq1 = 'ACTGTC'
seq2 = 'CTAG'



matrix = [[0 for x in range(len(seq1)+1)] for x in range(len(seq2)+1)]

# dynamically change width of each cell displayed
def add_space(string, width):
	
	string_len = len(string)
	while string_len < width:
		string += ' '
		string_len += 1
		
	return string

# print results in matrix format
score_lengths = []
for row in matrix:
	for score in row:
		score_lengths.append(len(str(score)))

width = max(score_lengths)

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
	


	

	
	
	
	
	
	
	
	
