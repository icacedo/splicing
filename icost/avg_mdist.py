import argparse
import json

parser = argparse.ArgumentParser()
parser.add_argument('json_file', type=str, metavar='<file>',
	help='input json file with tested icosts')

args = parser.parse_args()

f = open(args.json_file)

mdists = json.load(f)

for icost in mdists:
	print(icost[0])
	
