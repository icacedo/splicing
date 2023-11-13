import argparse
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('program', type=str, metavar='<program>',
	help='location of mdist program')
parser.add_argument('apc_dir', type=str, metavar='<apc_dir>',
	help='directory with apc generated gff files')
parser.add_argument('wb_dir', type=str, metavar='<apc_dir>',
	help='directory with wb gff files')

args = parser.parse_args()

apc_files = {}
for f in os.listdir(args.apc_dir):
	apath = args.apc_dir + f
	aID = f.split('.')[1]
	apc_files[aID] = apath

wb_files = {}
for f in os.listdir(args.wb_dir):
	wpath = args.wb_dir + f
	wID = f.split('.')[1]
	wb_files[wID] = wpath

IDcmds = {}
for ID in apc_files:
	prog = args.program
	agff = apc_files[ID]
	wgff = wb_files[ID]
	IDcmds[ID] = f'python3 {prog} {agff} {wgff}'.split(' ')

for ID in IDcmds:
	scmd = IDcmds[ID]
	proc = subprocess.Popen(scmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
		stdin=subprocess.PIPE, text=True)
	stdout, stderr = proc.communicate()
	print('Gene ID:', ID)
	print(stdout)
