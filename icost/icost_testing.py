import subprocess
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('apc_pkls', type=str, metavar='<directory>', 
	help='input directory with apc pickle files')
parser.add_argument('apc_fastas', type=str, metavar='<directory>',
	help='input directory with apc fasta files')

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
parser.add_argument('--acceptor_pwm', required=False, type=str, metavar='<file>',
	help='acceptor pwm .tsv')

parser.add_argument('--icost_range', required=False, type=float, default=100.00,
	metavar='<float>', help='intron cost %(default).2f')

args = parser.parse_args()

program = 'apc_score.py'
apc_dir = args.apc_fastas
pkl_dir = args.apc_pkls
outdir = 'icost_testing_out/'
os.makedirs(os.path.dirname(outdir), exist_ok=True)

exon_mm = args.exon_mm
intron_mm = args.intron_mm
exon_len = args.exon_len
intron_len = args.intron_len
donor_pwm = args.donor_pwm
acceptor_pwm = args.acceptor_pwm

fasta_paths = {}
for fname in os.listdir(apc_dir):
	if fname.endswith('fa'):
		ID1 = fname.split('.')[1]
		fpath = apc_dir + fname
		fasta_paths[ID1] = fpath

print(fasta_paths)

pkl_paths = {}
for fname in os.listdir(pkl_dir):
	ID2 = fname.split('.')[1]
	ppath = pkl_dir + fname
	pkl_paths[ID2] = ppath

print(pkl_paths)

'''
subprocess.run(f'python3 {apc_score} {pkl_file} {fa_file} {exon_mm}'
	f' {intron_mm} {exon_len} {intron_len} {donor_pwm} {acceptor_pwm}'
	f' > {outdir}', shell=True)

gff1 = 'outdir/test'
gff2 = 'data/build/apc/
'''


