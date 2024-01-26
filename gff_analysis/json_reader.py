import argparse
import os
import json

parser = argparse.ArgumentParser()
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
print('#####')
# next check if other isoforms matche wb

wbmatch2 = {}
wbmatch3 = {}
for fname in os.listdir(args.json_dir):
	gID = fname.split('.')[0]
	with open(args.json_dir+fname) as jfile:
		info = json.load(jfile)
		matches = []
		for iso in info:
			if iso == f'ch.{gID}-1' or iso == f'ch.{gID}-wb': continue
			if info[iso]['wb_frame'] == True:
				wbmatch = (iso, info[iso]['prob'])	
				matches.append(wbmatch)
			wbmatch3[gID] = wbmatch
		wbmatch2[gID] = matches

for gID in wbmatch2:
	print(gID, wbmatch2[gID])

print(len(wbmatch2))
print(len(wbmatch3))

for gID in wbmatch3:
	print(gID, wbmatch3[gID])
