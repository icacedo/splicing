import argparse
import isosort_lib as isl
import os

fasta = sys.argv[1]
jfile = sys.argv[2]

parser = argparse.ArgumentParser()
parser.add_argument('fa_dir', type=str, metavar='<directory>',
	help='directory with wormbase fasta files')
parser.add_argument('j_dir', type=str, metavar='<directory>',
	help='directory with apc isoform json files')
parser.add_argument('g_list', type=str, metavar='<file>',
	help='text file with list of bad isos to examine')

args = parser.parse_args()

seq = isl.get_seq(fasta)
seq_mod = isl.mod_seq(seq)
sym_seq_wb = isl.make_wb_sym(jfile, seq_mod)
sym_seq_apc = isl.make_apc_sym(jfile, seq_mod)
frame_sym = isl.make_frame_sym(sym_seq_wb)
	
print(frame_sym)
print(sym_seq_wb)
print(seq_mod)
print(sym_seq_apc)

