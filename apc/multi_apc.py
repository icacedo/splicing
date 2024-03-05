# based on Ian's code in setup/bin/parallelize
# caution: multi.py will run even if dependencies are not available
# output will be incorrect

import argparse
import multiprocessing as mp
import subprocess
import os
import time

parser = argparse.ArgumentParser(description='command line parallelizer')
parser.add_argument('file', help='file of command lines')
parser.add_argument('--cpus', required=False, type=int, default=1,
	metavar='<int>', help='number of CPUs to use [%(default)i]')

args = parser.parse_args()

def worker(cmd):
	fpath = cmd.split(' ')[2]
	gID = fpath.split('.')[1]
	print(f'working on {gID}')
	return subprocess.run(cmd, shell=True, capture_output=True).stdout.decode()

jobs = []
dpath = ''
with open(args.file, 'r') as fp:
	for cmd in fp:
		path = cmd.split('>')[1]
		path = path.split('/')[1:-1]
		dpath = '/' + '/'.join(path) + '/'	
		jobs.append(cmd.rstrip())

os.makedirs(dpath, exist_ok=True)

pool = mp.Pool(args.cpus)
pool.map(worker, jobs)

# testing time
'''	
starttime = time.time()
pool = mp.Pool(args.cpus)
for result in pool.map(worker, jobs):
	print(result, end='')
endtime = time.time()
print(endtime-starttime)

starttime = time.time()
for i in jobs:
	worker(i)
endtime = time.time()
print(endtime-starttime)
'''







