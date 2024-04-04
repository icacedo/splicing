import argparse
import csv
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('gene_list', help='csv file with list of genes to use')
parser.add_argument('apcgen_dir', help='directory with apcgen gff files')
parser.add_argument('--out_dir', help='name of output directory')
parser.add_argument('--iso_type', help=f'select which type of gene you want' 
					f' to create jbrowse gffs for; i.e. input 2 for genes'
					f' where the 2nd isoform matches wormbase; input 0 for'
					f' genes where no isoforms match wormbase')

args = parser.parse_args()

glist = []
with open(args.gene_list, newline='') as csvfile:
	reader = csv.DictReader(csvfile)
	for row in reader:
		if row['wb match'] == args.iso_type:
			glist.append(row['iso id'].split('.')[1])

fpaths = {}
for file in os.listdir(args.apcgen_dir):
	gid = file.split('.')[1]
	if gid in glist:
		fpaths[gid] = f'{args.apcgen_dir}{file}'

program1 = 'apcgen2jb.py'
program2 = 'splitgff.py'

if os.path.exists(args.out_dir): pass
else: os.mkdir(args.out_dir)

for gid in fpaths:	
	path = fpaths[gid]
	outf = args.out_dir + gid + '.jb.gff'
	subprocess.run(f'python3 {program1} {path} > {outf}', 
					shell=True)
	subprocess.run(f'python3 {program2} {outf} {args.out_dir}', 
					shell=True)
	print(f'{gid}.jb.gff')






