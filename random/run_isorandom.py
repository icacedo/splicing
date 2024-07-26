'''
How to setup up genomikon/isoformer/isorandom.c
Follow instructions on genomikon github to compile binaries
In Code/bin/ add symlink 
ln -s PATH/TO/genomikon/isoformer/isorandom
Add export PATH=$PATH:$HOME/Code/bin to .bashrc
Can now execute isorandom from anywhere
'''

import argparse
import subprocess
import multiprocessing 

parser = argparse.ArgumentParser('wrapper for isorandom.c')

parser.add_argument('beg_len', type=int, metavar='<int>',
    help='length of shortest sequence to test')
parser.add_argument('end_len', type=int, metavar='<int>',
    help='length of longest sequence to test')
parser.add_argument('inc', type=int, metavar='<int>',
    help='amount to increment tested sequence lengths by')
parser.add_argument('count', type=int, metavar='<int>',
    help='number of simulations to run')
parser.add_argument('--cpus', required=False, type=int, default=1,
    metavar='<int>', help='number of CPUs to use [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25, 
    metavar='<int>', help='minimum exon length [%(default)i]')
parser.add_argument('--min_intron', required=False, type=int, default=35, 
    metavar='<int>', help='minimum intron length [%(default)i]')
parser.add_argument('--max_splice', required=False, type=int, default=3, 
    metavar='<int>', help='maximum splices [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=100, 
    metavar='<int>', help='genomic flank lengths [%(default)i]')
# --seed does not work as an arg for isorandom
parser.add_argument('--seed', required=False, type=int, default=1,
    metavar='<int>', help='set random seed [%(default)i]')

args = parser.parse_args()

def worker(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, 
            text=True)
'''
with open('isorandom.out', 'w') as file:
    for i in range(args.beg_len, args.end_len + args.inc, args.inc):
        result = subprocess.run([
            "isorandom", str(i), str(args.count), "--min_exon", str(args.min_exon),
            "--min_intron", str(args.min_intron), "--max_splice", 
            str(args.max_splice), "--flank", str(args.flank)
        ], check=True, text=True, capture_output=True)
        file.write(result.stdout)
'''

jobs = []
for i in range(args.beg_len, args.end_len + args.inc, args.inc):
    cmd = (
        f"isorandom {i} {args.count} --min_exon {args.min_exon}"
        f"--min_intron {args.min_intron} --max_splice {args.max_splice}"
        f"--flank {args.flank} > {i}.out"
    )
    jobs.append(cmd)

pool = multiprocessing.Pool(args.cpus)
pool.map(worker, jobs)


