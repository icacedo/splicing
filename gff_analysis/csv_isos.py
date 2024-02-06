import argparse
import os
import json
import csv

parser = argparse.ArgumentParser()
parser.add_argument('json_dir', type=str, metavar='<directory>',
	help='directory with isosort json files')

args = parser.parse_args()

# only write top isos
with open('apcgen_isos.csv', 'w', newline='') as csvfile:
	fieldnames = ['gID', 'WBGene']
	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
	writer.writeheader()

for fname in os.listdir(args.json_dir):
	gID = fname.split('.')[0]
	with open(args.json_dir+fname) as jfile:
		info = json.load(jfile)
		for iso in info:
			if iso == f'ch.{gID}-1' and info[iso]['wb_frame'] == True:
				print(iso, info[iso])
	



