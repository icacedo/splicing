import sys
import csv

cfile = sys.argv[1]

genes = []
with open(cfile, newline='') as csvfile:
	reader = csv.DictReader(csvfile)
	for row in reader:
		genes.append(row)
		
by_fit = sorted(genes, key=lambda d: float(d['fitness']))

print(f'gid,fitness,wbgene,wbmatch')

# get best predictions
'''
for f in by_fit:
	if f['wbmatch'] == '1':
		print(f'{f['gid']},{f['fitness']},{f['wbgene']},{f['wbmatch']}')


# get 2nd best predictions

for f in by_fit:
	if f['wbmatch'] == '2':
		print(f'{f['gid']},{f['fitness']},{f['wbgene']},{f['wbmatch']}')


# get bad predictions
for f in by_fit:
	if f['wbmatch'] == 'None': continue
	if int(f['wbmatch']) <= 10:
		print(f'{f['gid']},{f['fitness']},{f['wbgene']},{f['wbmatch']}')
'''	
# get worst predictions
for f in by_fit:
	if f['wbmatch'] == 'None': 
		print(f'{f['gid']},{f['fitness']},{f['wbgene']},{f['wbmatch']}')

