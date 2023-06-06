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
parser.add_argument('-mm', action='store_true')
parser.add_argument('-len', action='store_true')
parser.add_argument('-pwm', action='store_true')

args = parser.parse_args()

#if args.mm:
#	print('wowowo')

#exins = ml.read_txt_seqs(args.extxt)

def mm_tsv_write(exins):

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

if args.extxt:
	exons = ml.read_txt_seqs(args.extxt)
if args.intxt:
	introns = ml.read_txt_seqs(args.intxt)

# this works but it might be too confusing
if args.extxt and args.mm:
	print('mm only ex')
elif args.extxt and args.len:
	print('len only ex')
elif args.extxt:
	print('mm and len ex')

if args.intxt and args.mm:
	print('mm only in')
elif args.intxt and args.len:
	print('len only in')
elif args.intxt:
	print('mm and len in')


'''
if args.extxt and args.mm:
	mm_tsv_write(exins)
if args.intxt and args.mm:
	mm_tsv_write(exins)
'''











