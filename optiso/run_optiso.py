import argparse
import json
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('config_dir', type=str, 
	help='directory with config.json files for individual genes')
parser.add_argument('--cpu', required=False, type=int, default=1, 
	help='number of cpus to use')

arg = parser.parse_args()

sum_params = {}
for file in os.listdir(arg.config_dir):
	name = f'ch.{file.split('.')[0]}'
	cmd = f'./optiso {arg.config_dir}{file} --cpu {arg.cpu}'
	result = subprocess.run(cmd, shell=True, capture_output=True)
	jstring = result.stdout.decode('utf-8')
	ginfo = json.loads(jstring)
	sum_params[name] = ginfo

print(sum_params)

	
