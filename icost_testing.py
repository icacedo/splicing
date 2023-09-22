import subprocess
import os

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




