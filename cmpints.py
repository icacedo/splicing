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

cmds = {}
for ID in apc_files:
	cmd = args.program
	agff = apc_files[ID]
	wgff = wb_files[ID]
	print(f'python3 {cmd} {agff} {wgff}')
	






