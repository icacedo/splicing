import isosort_lib as isl
import argparse
import os
import json

parser = argparse.ArgumentParser(
	description='adds labels to apc generated isoforms,' 
				' writes to .json files in out/')
parser.add_argument('wb_dir', type=str, metavar='<directory>',
	help='directory with wormbase fasta/gff files')
parser.add_argument('apcgen_dir', type=str, metavar='<directory>',
	help='directory with apc generated gff files')
parser.add_argument('--weights', type=str, metavar='<file>',
	help='file with individual gene weights')
parser.add_argument('--elen', type=str, metavar='<file>',
	help='exon length model .tsv')
parser.add_argument('--ilen', type=str, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--emm', type=str, metavar='<file>',
	help='exon length model .tsv')
parser.add_argument('--imm', type=str, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--dpwm', type=str, metavar='<file>',
	help='donor site pwm .tsv')
parser.add_argument('--apwm', type=str, metavar='<file>',
	help='acceptor site pwm .tsv')

args = parser.parse_args()

fastas = []
gffs = []
for fname in os.listdir(args.wb_dir):
	if fname.endswith('fa'): fastas.append(fname)
	if fname.endswith('gff3'): gffs.append(fname)

paths = {}
for fname in fastas:
	gID = fname.split('.')[1]
	paths[gID] = [args.wb_dir+fname]

for fname in gffs:
	gID = fname.split('.')[1]
	paths[gID] += [args.wb_dir+fname]

for fname in os.listdir(args.apcgen_dir):
	gID = fname.split('.')[1]
	paths[gID] += [args.apcgen_dir+fname]

gwts = {}
with open(args.weights, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		line = line.split('\t')
		fit = line[0]
		wdpwm = line[1]
		wapwm = line[2]
		wemm = line[3]
		wimm = line[4]
		welen = line[5]
		wilen = line[6]
		icost = line[7]
		gid = line[8]
		gwts[f'ch.{gid}'] = {
			'fit': fit,
			'wdpwm': wdpwm,
			'wapwm': wapwm,
			'wemm': wemm,
			'wimm': wimm,
			'welen': welen,
			'wilen': wilen
			}

os.makedirs('sort_out/', exist_ok=True)

for gID in paths:
	if gID == 'ch.4567': continue
	fa = paths[gID][0]
	wb_gff = paths[gID][1]
	abcgen_gff = paths[gID][2]
	seq = isl.get_seq(fa)
	wbginfo = isl.get_wbgene_info(wb_gff, seq)
	wbginfo = isl.score_wb_iso(
		seq, wbginfo, 
		args.elen, args.ilen, args.emm, 
		args.imm, args.dpwm, args.apwm, 
		gwts[gID]['welen'], gwts[gID]['wilen'], gwts[gID]['wemm'], 
		gwts[gID]['wimm'], gwts[gID]['wdpwm'], gwts[gID]['wapwm'],
		gwts[gID]['icost']
	)
	print(wbginfo)
	break


'''
for gID in paths:
	fa = paths[gID][0]
	wb_gff = paths[gID][1]
	apcgen_gff = paths[gID][2]
	isoforms_info = isl.amass_info(fa, wb_gff, apcgen_gff, args.elen, 
									args.ilen, args.emm, args.imm, 
									args.dpwm, args.apwm, gwts[gid])
	jstr = json.dumps(isoforms_info, indent = 4)
	with open(f'out/{gID}.json', 'w') as outfile:
		outfile.write(jstr)
'''









