import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('json_file', type=str, metavar='<file>',
	help='input json file with tested icosts')

args = parser.parse_args()

f = open(args.json_file)

mdists = json.load(f)

avg_mdists = {}
for icost in mdists:
	avg_mdists[icost[0]] = 0
	for mdist_dict in icost[1]:
		avg_mdists[icost[0]] += float(mdist_dict['mdist'])/len(icost[1])

for icost in avg_mdists:
	print('icost:', icost, 'avg mdist:', avg_mdists[icost])
