import argparse
import csv
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('gene_list', help='csv file with list of genes to use')
parser.add_argument('apcgen_dir', help='directory with apcgen gff files')
parser.add_argument('--iso_type', help=f'select which type of gene you want' 
					f' to create jbrowse gffs for; i.e. input 2 for genes'
					f' where the 2nd isoform matches wormbase; input 0 for'
					f' genes where no isoforms match wormbase')

args = parser.parse_args()

#redict = {}
with open(args.gene_list, newline='') as csvfile:
	reader = csv.DictReader(csvfile)
	for row in reader:
		if row['wb match'] == args.iso_type:
			print(row['iso id'], row['wb match'])
		#redict[row['iso id']] = row['wb match']

os.listdir(args.apcgen_dir)




