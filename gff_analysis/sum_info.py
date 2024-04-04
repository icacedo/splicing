import argparse
import os
import json
import csv

parser = argparse.ArgumentParser(
		description='generates lists of isoforms based on categories')
parser.add_argument('json_dir', type=str, metavar='<directory>',
	help='directory with isosort json files')

args = parser.parse_args()

total = len(os.listdir(args.json_dir))
wbmatch = []
for fname in os.listdir(args.json_dir):
	gID = fname.split('.')[0]
	with open(args.json_dir+fname) as jfile:
		info = json.load(jfile)
		if info[f'ch.{gID}-1']['wb_frame'] == True:
			prob = info[f'ch.{gID}-1']['prob']
			wbmatch.append((fname, prob))

wbmatch_bins = {
	0.99: [],
	0.90: [],
	0.80: [],
	0.70: [],
	0.60: [],
	0.50: [],
	0.40: [],
	0.30: [],
	0.20: [],
	0.10: [],
	0.00: []
}
for tup in wbmatch:
	if tup[1] >= 0.99:
		wbmatch_bins[0.99].append(tup)
		continue
	if tup[1] >= 0.9 and tup[1] < 0.99:
		wbmatch_bins[0.9].append(tup)
		continue
	for i in reversed(range(0, 9, 1)):
		if tup[1] >= i/10 and tup[1] < round((i/10)+0.1, 1):
			wbmatch_bins[round(i*0.1, 1)].append(tup)

print(f'top iso match: {len(wbmatch)} out of {total}')
print(f'bin\tiso count')
for bn in wbmatch_bins:
	print(f'{bn}\t{len(wbmatch_bins[bn])}')

gmatches = []
for gene in wbmatch:
	gmatches.append(gene[0].split('.')[0])

wbmatch2 = []
for fname in os.listdir(args.json_dir):
	gID = fname.split('.')[0]
	if gID in gmatches: continue
	with open(args.json_dir+fname) as jfile:
		info = json.load(jfile)
	for iso in info:
		if iso == f'ch.{gID}-1' or iso == f'ch.{gID}-wb': continue
		if info[iso]['wb_frame'] == True:
			match2 = (iso, info[iso]['prob'])
			wbmatch2.append(match2)		

gmatches2 = []
for gene in wbmatch2:
	iid = gene[0].split('.')[1]
	gmatches2.append(iid.split('-')[0])

nomatch = []
for fname in os.listdir(args.json_dir):
	gID = fname.split('.')[0]
	if gID in gmatches: continue
	if gID in gmatches2: continue
	nomatch.append(gID)

total_nomatch = len(wbmatch2) + len(nomatch)
print(f'top iso no match: {total_nomatch} out of {total}')
print(f'other iso match: {len(wbmatch2)} out of {total_nomatch}' )

match2nd = []
for w in wbmatch2:
	if w[0].split('-')[1] == str(2):
		match2nd.append(w[0].split('-')[0])

print(f'2nd iso match: {len(match2nd)} out of {len(wbmatch2)}')
print(f'no iso match: {len(nomatch)} out of {total_nomatch}')

isomatches = {}
iso_wbgenes = {}
for fname in os.listdir(args.json_dir):
	gID = fname.split('.')[0]
	with open(args.json_dir+fname) as jfile:
		info = json.load(jfile)
	has_match = False
	for iso in info:
		iid = iso.split('-')
		if 'wb' in iso: 
			iso_wbgenes[iid[0]] = info[iso]['Parent=Gene']
			continue
		if info[iso]['wb_frame'] == True:
			has_match = True
			isomatches[iid[0]] = int(iid[1])
	if has_match == False:
		isomatches[f'ch.{gID}'] = 0

dlist = []
for iid in isomatches:
	d = {}
	d['iso id'] = iid
	d['wb match'] = isomatches[iid]
	d['WBGene'] = iso_wbgenes[iid]
	dlist.append(d)

fields = ['iso id', 'WBGene', 'wb match']
filename = 'wb_matching_isos.csv'
with open(filename, 'w') as csvfile:
	writer = csv.DictWriter(csvfile, fieldnames=fields)
	writer.writeheader()
	writer.writerows(dlist)

		

