import sys
import csv

# 1045info.csv
cfile = sys.argv[1]

genes = []
with open(cfile, newline='') as csvfile:
	reader = csv.DictReader(csvfile)
	for row in reader:
		genes.append(row)

'''
match_counts = {}
for f in genes:
    m = f['wbmatch']
    if m in match_counts:
        match_counts[m] += 1
    else:
        match_counts[m] = 0
        match_counts[m] += 1
        
print(f'wbmatch,counts')
for m in match_counts:
    if m == 'None':
        print(f'0.0,{match_counts[m]}')
    else:
        print(f'{m},{match_counts[m]}')


        #print(f'{f['gid']},{f['fitness']},{f['wbgene']},{f['wbmatch']}')
'''
print('fit')
fits = []
for f in genes:
	fits.append(f['fitness'])
for f in fits:
	print(f)

