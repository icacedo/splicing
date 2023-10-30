# based on Ian's code in setup/bin/parallelize

import argparse
import multiprocessing as mp
import subprocess
import os

parser = argparse.ArgumentParser(description='command line parallelizer')
parser.add_argument('file', help='file of command lines')
parser.add_argument('--cpus', required=False, type=int, default=1,
	metavar='<int>', help='number of CPUs to use [%(default)i]')

args = parser.parse_args()


def worker(cmd):
	return subprocess.run(cmd, shell=True, capture_output=True).stdout.decode()

jobs = []
with open(args.file, 'r') as fp:
	for cmd in fp:
		jobs.append(cmd.rstrip())

# this puts results in a list
# can i do that for the resulting gffs, and put them in individual files?
pool = mp.Pool(args.cpus)
a = pool.map(worker, jobs)
print(a)
for i in a:
	print(i)

	
'''
with Pool(5) as p:
	print(p.map(f, [1, 2, 3]))

# go to codesdope.com intro to multiprocessing and process in python

print(os.cpu_count())
'''





