import argparse
import os
import subprocess

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
	wb_gffs[wID] = wb_path

def find_mdist(cap):
		
	dcap = cap.stdout.decode('UTF-8')
	for line in dcap.split('\n'):
		if line.startswith('mdist'):
			mdist = line.split('=')[1]
			return mdist

mdist_groups = {}
for ID in icost_gffs:
	for path in icost_gffs[ID]:
		gff1_apc = path 
		gff2_wb = wb_gffs[ID]
		icost = gff1_apc.split('_')[2]
		print(ID)
		print(icost)
		print(gff1_apc)
		print(gff2_wb)
		cap = subprocess.run(f'python3 {program} {gff1_apc} {gff2_wb}', 
			shell=True, capture_output=True)
		print(find_mdist(cap))
		break
	break
