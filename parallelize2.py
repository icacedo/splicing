#!/usr/bin/env python3

# based on Ian's code in setup/bin/parallelize

import argparse
from multiprocessing import Pool
import subprocess
import os

def worker(cmd):
	return subprocess.run(cmd, shell=True, capture_output=True).stdout.decode()

def f(x):
	return x*x

with Pool(5) as p:
	print(p.map(f, [1, 2, 3]))

from multiprocessing import Process

# go to codesdope.com intro to multiprocessing and process in python

print(os.cpu_count())






