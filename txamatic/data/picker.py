import os
import sys
import json

fields = ['region', 'length', 'introns', 'strand', 'RNASeq_splice', 'max_splices']
print('\t'.join(fields))

for id in os.listdir(sys.argv[1]):
	with open(f'{sys.argv[1]}/{id}/{id}.json') as fp:
		s = fp.read()
		j = json.loads(s)
		
		if j['strand'] != '+': continue
		if j['length'] > 1000: continue
		if j['introns'] < 2: continue
		if j['max_splices'] < 100000: continue
		
		for k in fields:
			print(j[k], end="\t")
		print()

