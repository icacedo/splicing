import argparse
import gzip
import modelib as ml

parser = argparse.ArgumentParser(
	description='Generates len, MM, and PWM models for apc')

parser.add_argument('--extxt', type=str, metavar='<file>', 
	help='input text file with exon sequences')
parser.add_argument('--intxt', type=str, metavar='<file>',
	help='input text file with intron sequences')
parser.add_argument('--dntxt', type=str, metavar='<file>',
	help='input text file with donor site sequences')
parser.add_argument('--actxt', type=str, metavar='<file>',
	help='input text file with acceptor site sequences')
args = parser.parse_args()

# length model

exons = ml.read_txt_seqs(args.extxt)

exlen_scores, exlen_probs = ml.memoize_fdist(exons)
print(exlen_scores)

exmm_scores, exmm_probs = ml.make_mm(exons)
print(exmm_scores)
print('**********')
print(exmm_probs)



