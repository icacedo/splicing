import argparse
import random
import json

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
		help='input single gene fasta file')
parser.add_argument('gff3', type=str, metavar='<file>',
		help='input single gene gff3 file')
parser.add_argument('--program', require=True, type=str, 
		metavar='<exec>'

parser.add_ar





parser.add_argument('config', type=str, metavar='<json>',
		help='configuration file')

parser.add_argument('--pop', required=False, type=int, default=50,
		metavar='<int>', help='population size [%(default)i]')
parser.add_argument('--seed', required=False, type=int,
		metavar='<int>', help='random seed')

args = parser.parse_args()

def read_cfg(filename):
	
		cfg = None
		with open(filename) as fp:
			s = fp.read()
			cfg = json.loads(s)

		return cfg


CFG = read_cfg(args.config)

print(CFG)






if args.seed: random.seed(args.seed)







