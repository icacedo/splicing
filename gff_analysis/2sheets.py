import sys
import csv

# write list of genes into csv for google sheets
file = sys.argv[1]

genes = []
with open(file, 'r') as csvfile:
	reader = csv.DictReader(csvfile)
	d = {}
	for row in reader:
		if row['wb match'] == '2':
			d = {}
			d['iso id'] = row['iso id']
			d['WBGene'] = row['WBGene']
			d['wb match'] = row['wb match']
			genes.append(d)

fields = ['iso id', 'WBGene', 'wb match']
with open('2isos.csv', 'w') as csvfile:
	writer = csv.DictWriter(csvfile, fieldnames=fields)
	writer.writeheader()
	writer.writerows(genes)
				
