import subprocess
import os
import argparse
import apc_score.py as apcsc

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

'''
apc_score = 'apc_score.py'
fa_file = 'data/build/apc/ch.9940.fa'
pkl_file = 'apc_isoforms.pkl'
outdir = 'icost_testing_out/'
os.makedirs(os.path.dirname(outdir), exist_ok=True)

exon_mm = '--exon_mm mkmdls_out/exon_mm.tsv'
intron_mm = '--intron_mm mkmdls_out/intron_mm.tsv'
exon_len = '--exon_len mkmdls_out/exon_len.tsv'
intron_len = '--intron_len mkmdls_out/intron_len.tsv'
donor_pwm = '--donor_pwm mkmdls_out/donor_pwm.tsv'
acceptor_pwm = '--acceptor_pwm mkmdls_out/acceptor_pwm.tsv'

subprocess.run(f'python3 {apc_score} {pkl_file} {fa_file} {exon_mm}'
	f' {intron_mm} {exon_len} {intron_len} {donor_pwm} {acceptor_pwm}'
	f' > {outdir}test', shell=True
)

gff1 = 'outdir/test'
gff2 = 'data/build/apc/
'''


