import argparse
import os
import csv
import json

parser = argparse.ArgumentParser()
parser.add_argument('glist', type=str, metavar='<file>',
	help='txt file with genes to use')
parser.add_argument('jdir', type=str, metavar='<directory>',
	help='directory with json files')

args = parser.parse_args()
'''
with open('match2_isos.csv', 'w', newline='') as csvfile:
	fieldnames = ['gID', 'WBGene', 'top_iso_prob', 
'''
gIDs = []
with open(args.glist, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		gIDs.append(line)

badiso_info = {}
for j in os.listdir(args.jdir):
	jID = 'ch.' + j.split('.')[0]
	if jID in gIDs:
		with open(args.jdir+j) as jfile:
			info = json.load(jfile)	
			badisos = []
			for gID in info:
				if gID.endswith('-1') or gID.endswith('-wb'):
					badisos.append(info[gID])
			badiso_info[jID] = badisos

with open('badisos2.csv', 'w', newline='') as csvfile:
	cw = csv.writer(csvfile, delimiter=',', quotechar='"')
	cw.writerow(['ID', 'WBGene', 'prob', 'exons', 'introns', 
		'in_frame', 'codons', 'PTC'])	
	for gID in badiso_info:
		apc = badiso_info[gID][0]
		wb = badiso_info[gID][1]
		wbg = wb['Parent=Gene']
		prob = apc['prob']
		exons = apc['exons']
		introns = apc['introns']
		in_frame = apc['in_frame']
		codons = apc['codons']
		PTC = apc['PTC']
		cw.writerow([gID, wbg, prob, exons, introns, in_frame, codons,
			PTC])
		print(apc)
		print('#####')
		print(wb)
		break
