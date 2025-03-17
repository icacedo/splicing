import argparse

parser = argparse.ArgumentParser()

parser.add_argument('apc_isos', type=str, metavar='<directory>',
	help='directory with APC generated gff files')

args = parser.parse_args()