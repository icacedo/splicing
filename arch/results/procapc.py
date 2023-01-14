import json
import sys


dset = {}

with open(sys.argv[1]) as fp:
	for line in fp.readlines():
		if len(line) < 1: continue
		name, js = line.split('\t')
		data = json.loads(js)
		dset[name] = {
			'don': data['genotype']['--wdpwm'],
			'acc': data['genotype']['--wapwm'],
			'emm': data['genotype']['--wemm'],
			'imm': data['genotype']['--wimm'],
			'elen': data['genotype']['--welen'],
			'ilen': data['genotype']['--wilen'],
			'cost': data['genotype']['--icost'],
			'fit': data['fitness'],
		}

fields = ['don', 'acc', 'emm', 'imm', 'elen', 'ilen', 'cost', 'fit']
'''
for field in fields:
	minv = 1e9
	maxv = -1e9
	for ch in dset:
		if dset[ch][field] < minv: minv = dset[ch][field]
		if dset[ch][field] > maxv: maxv = dset[ch][field]
	sumv = maxv - minv
	for ch in dset:
		dset[ch][field] = int(100 * (dset[ch][field] - minv) / sumv)
'''

print('seq', end='\t')
print('\t'.join(fields))
for ch in dset:
	print(ch, end='')
	for field in fields:
		print(f'\t{dset[ch][field]}', end='')
	print()
