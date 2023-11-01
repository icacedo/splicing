# based on Ian's code in setup/bin/parallelize

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
	return subprocess.run(cmd, shell=True, capture_output=True).stdout.decode()

jobs = []
with open(args.file, 'r') as fp:
	for cmd in fp:
		jobs.append(cmd.rstrip())

def square(x):
	print(f'start process:{x}')
	square = x * x
	print(f'square {x}:{square}')
	time.sleep(1)
	print(f'end process{x}')

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
starttime = time.time()
pool = mp.Pool()
pool.map(square, range(0, 5))
pool.close
endtime = time.time()
print(f'Time taken {endtime-starttime} seconds')

print('**********')

def s2(x):
	square = x * x
	print('start process')
	time.sleep(1)
	return square

pool = mp.Pool()
result = pool.map_async(s2, range(0, 5))
print('main script')
print(result.get())
print('end main script')

pool = mp.Pool()
result = pool.map(s2, range(0, 5))
for i in result:
	print(i)
'''





'''
nums = []
for i in range(0, 5):
	nums.append(i)
print(nums)

pool = mp.Pool(args.cpus)
a = pool.map(jobs, nums)
for i in a:
	print(i)
'''
'''
pool = mp.Pool(args.cpus)

for result in pool.map(worker, jobs):
	print(result, end='')
'''

'''
# this puts results in a list
# can i do that for the resulting gffs, and put them in individual files?
pool = mp.Pool(args.cpus)
a = pool.map(worker, jobs)
print(a)
for i in a:
	print(i)
'''
	
'''
with Pool(5) as p:
	print(p.map(f, [1, 2, 3]))

# go to codesdope.com intro to multiprocessing and process in python

print(os.cpu_count())
'''





