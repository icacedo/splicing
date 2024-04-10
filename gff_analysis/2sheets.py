import sys
import csv

# write list of genes into csv for google sheets
file = sys.argv[1]

redict = {}
with open(file, 'r') as csvfile:
	reader = csv.DictReader(csvfile)
	for row in reader:
		ginfo = {}
		ginfo['iso id'] = row['iso id']
		ginfo['WBGene'] = row['WBGene']
		if int(row['wb match']) not in redict:
			redict[int(row['wb match'])] = []
			redict[int(row['wb match'])].append(ginfo)
		else:
			redict[int(row['wb match'])].append(ginfo)

redict = dict(sorted(redict.items()))
fields = ['iso id', 'WBGene', 'wb match'] 
for n in redict:
	with open(f'{n}isos.csv', 'w') as csvfile:
		genes = []
		for g in redict[n]:
			d = {}	
			d['iso id'] = g['iso id'] 
			d['WBGene'] = g['WBGene']
			d['wb match'] = n
			genes.append(d)		
		writer = csv.DictWriter(csvfile, fieldnames=fields)
		writer.writeheader()
		writer.writerows(genes)
		csvfile.close()
		
		












