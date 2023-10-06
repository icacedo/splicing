import subprocess
import os
import argparse
import numpy as np
import mdist_lib as mdl

parser = argparse.ArgumentParser()
parser.add_argument('apc_pkls', type=str, metavar='<directory>', 
	help='input directory with apc pickle files')
parser.add_argument('apc_dir', type=str, metavar='<directory>',
	help='input directory with apc fasta gff files')
parser.add_argument('--path2ml', type=str, metavar='<directory path>', 
	help='absolute path to directory with modelib')
parser.add_argument('--outdir', type=str, metavar='<outdir path>',
	required=True, help='/path/ to outdir')

parser.add_argument('--exon_len', type=str, metavar='<file>', 
	help='exon length model .tsv')
parser.add_argument('--intron_len', type=str, metavar='<file>',
	help='intron length model .tsv')
parser.add_argument('--exon_mm', type=str, metavar='<file>',
	help='exon markov model .tsv')
parser.add_argument('--intron_mm', type=str, metavar='<file>',
	help='intron markov model .tsv')
parser.add_argument('--donor_pwm', type=str, metavar='<file>',
	help='donor pwm .tsv')
parser.add_argument('--acceptor_pwm', type=str, metavar='<file>',
	help='acceptor pwm .tsv')

parser.add_argument('--icost_range', required=False, type=float, default=100,
	metavar='<int>', help='intron cost %(default)i')
parser.add_argument('--icost_step', required=False, type=float, default=0.1,
	metavar='<float>', help='intron cost step %(default).1f')

args = parser.parse_args()

program = 'apc_score.py'
apc_dir = args.apc_dir
pkl_dir = args.apc_pkls
outdir = args.outdir+'icost_out/'
os.makedirs(os.path.dirname(outdir), exist_ok=True)

exon_mm = args.exon_mm
intron_mm = args.intron_mm
exon_len = args.exon_len
intron_len = args.intron_len
donor_pwm = args.donor_pwm
acceptor_pwm = args.acceptor_pwm
mlpath = args.path2ml

fasta_paths = {}
for fname in os.listdir(apc_dir):
	if fname.endswith('fa'):
		ID1 = fname.split('.')[1]
		fpath = apc_dir + fname
		fasta_paths[ID1] = fpath

pkl_paths = {}
for fname in os.listdir(pkl_dir):
	ID2 = fname.split('.')[1]
	ppath = pkl_dir + fname
	pkl_paths[ID2] = ppath

irange = int(args.icost_range)
irange_step = args.icost_step

wb_gffs = {}
for wbfile in os.listdir(apc_dir):
	if wbfile.endswith('.fa'): continue
	wID = wbfile.split('.')[1]
	wb_path = apc_dir + wbfile
	wb_gffs[wID] = wb_path

for i in np.arange(0, irange+0.1, irange_step):
	icost = round(i, 1)
	for ID in pkl_paths:
		pkl_file = pkl_paths[ID]
		fa_file = fasta_paths[ID]
		if 'bli' in pkl_paths[ID].split('.'):
			gff_name = 'ch.'+ID+'.icost_'+str(icost)+'_'+'bli.gff'
		else:
			gff_name = 'ch.'+ID+'.icost_'+str(icost)+'_'+'apc.gff'
		print(gff_name)
		print(wb_gffs[ID])
		subprocess.run(f'python3 {program} {pkl_file} {fa_file}'
			f' --exon_len {exon_len} --intron_len {intron_len}'
			f' --exon_mm {exon_mm} --intron_mm {intron_mm}'
			f' --donor_pwm {donor_pwm} --acceptor_pwm {acceptor_pwm}'
			f' --path2ml {mlpath} --icost {icost} > {outdir}{gff_name}', 
				shell=True )
		print(outdir+gff_name)
		introns1 = mdl.get_gff_intron_probs(outdir+gff_name)
		introns2 = mdl.get_gff_intron_probs(wb_gffs[ID])
		mdist = mdl.get_mdist(introns1, introns2)
		print(mdist)
		print('#####')
