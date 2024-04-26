import argparse
import json
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('config_dir', type=str, 
	help='directory with config.json files for individual genes')
parser.add_argument('--program', required=False, type=str, 
	default='./isoformer', help='algorithm version to use %(default)s')
parser.add_argument('--cpu', required=False, type=int, default=1, 
	help='number of cpus to use %(default)d')

arg = parser.parse_args()

sum_params = {}
for file in os.listdir(arg.config_dir):
	iid = file.split('.')[0]
	name = f'ch.{iid}'	
	cmd = (
		f'./optiso {arg.config_dir}{file} --program ./isoformer '
		f'--cpu {arg.cpu}'
	)
	print(name)
	result = subprocess.run(cmd, shell=True, capture_output=True)
	print(result)
	result = subprocess.run(cmd, shell=True, capture_output=True)
	print(result)
	jstring = result.stdout.decode('utf-8')
	print(jstring)
	# fails make json object, idk why
	#ginfo = json.loads(jstring)
	#sum_params[name] = ginfo

#with open('optiso_out.json', 'w') as jfile:
#	jfile.write(json.dumps(sum_params, indent=4))
# fails 
	
