import argparse
import json
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('config_dir', type=str, 
	help='directory with config.json files for individual genes')
parser.add_argument('--program', required=False, type=str, 
	default='./isoformer', help='algorithm version to use %(default)s')

arg = parser.parse_args()

sum_params = {}
for file in os.listdir(arg.config_dir):
	iid = file.split('.')[0]
	name = f'ch.{iid}'	
	cmd = (
		f'./optiso {arg.config_dir}{file} --program {arg.program} '
	)
	print(name)
	result = subprocess.run(cmd, shell=True, capture_output=True)
	jstring = result.stdout.decode('utf-8')
	print(jstring)
	# fails make json object, idk why
	#ginfo = json.loads(jstring)
	#sum_params[name] = ginfo

#with open('optiso_out.json', 'w') as jfile:
#	jfile.write(json.dumps(sum_params, indent=4))
# fails 
	
