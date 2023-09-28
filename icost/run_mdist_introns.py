import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('apc_gffs', type=str, metavar='<directory>',
	help='input directory with apc generated gff files')
parser.add_argument('wb_gffs', type=str, metavar='<directory>',
	help='input directory with wb annotation gff files')

args = parser.parse_args()

program = 'mdist_introns.py'
apc_gff_dir = args.apc_gffs
wb_gff_dir = args.wb_gffs

icost_gffs = {}
for agff in os.listdir(apc_gff_dir):
	aID = agff.split('.')[1]	
	apc_path = apc_gff_dir + agff 
	if aID not in icost_gffs:
		icost_gffs[aID] = [apc_path]
	else:
		icost_gffs[aID] += [apc_path] 

wb_gffs = {}
for wgff in os.listdir(wb_gff_dir):
	if wgff.endswith('.fa'): continue
	wID = wgff.split('.')[1]
	wb_path = wb_gff_dir + wgff
	if wID not in icost_gffs: continue
	if wID not in wb_gffs:
		wb_gffs[wID] = [wb_path]
	else:
		wb_gffs[wID] += [wb_path]

print(icost_gffs)
print(wb_gffs)

