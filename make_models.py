import argparse
import gzip
import modelib as ml
import csv

parser = argparse.ArgumentParser(
	description='Generates len, MM, and PWM models for apc')

parser.add_argument('--extxt', type=str, metavar='<file>', 
	required=False, help='input text file with exon sequences')
parser.add_argument('--intxt', type=str, metavar='<file>',
	help='input text file with intron sequences')
parser.add_argument('--dntxt', type=str, metavar='<file>',
	help='input text file with donor site sequences')
parser.add_argument('--actxt', type=str, metavar='<file>',
	help='input text file with acceptor site sequences')
args = parser.parse_args()

def mm_tsv_write(exintxt):

	exins = ml.read_txt_seqs(exintxt)

	exinmm_scores, exinmm_probs = ml.make_mm(exins)

	root, ext = args.extxt.split('.')
	filename = root + '_mm' + '.tsv'

	with open(filename, 'w', newline='') as tsvfile:
		writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
		for key in exinmm_scores:
			writer.writerow([key + 'A', exinmm_probs[key][0], 
				exinmm_scores[key][0]])
			writer.writerow([key + 'C', exinmm_probs[key][1], 
				exinmm_scores[key][1]])
			writer.writerow([key + 'G', exinmm_probs[key][2], 
				exinmm_scores[key][2]])
			writer.writerow([key + 'T', exinmm_probs[key][3], 
				exinmm_scores[key][3]])
	tsvfile.close()

# input exon.txt output exon_mm.tsv
#mm_tsv_write(args.extxt)
# input intron.txt output intron_mm.tsv
#mm_tsv_write(args.intxt)

if args.extxt:
	mm_tsv_write(args.extxt)
if args.intxt:
	mm_tsv_write(args.intxt)
