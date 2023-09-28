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
listy = []
count = 0
for gff in os.listdir(apc_gff_dir):
	ID = gff.split('.')[1]	
	print(ID, type(ID))
	gpath = apc_gff_dir + gff
	if count == 0:
		icost_gffs[ID] = [gpath]
	else:
		icost_gffs[ID] += [gpath] 
	count += 1
	if count == 2: break
print(icost_gffs)


dish = {}
dish['ID'] = [1]
dish['ID'] += [2]
print(dish)




