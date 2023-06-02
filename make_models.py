import argparse
import gzip
import modelib as ml
import csv

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

exons = ml.read_txt_seqs(args.extxt)

#exlen_scores, exlen_probs = ml.memoize_fdist(exons)

exmm_scores, exmm_probs = ml.make_mm(exons)

for key in exmm_scores:
	print(key + 'A', exmm_probs[key][0], exmm_scores[key][0])
	print(key + 'C', exmm_probs[key][1], exmm_scores[key][1])
	print(key + 'G', exmm_probs[key][2], exmm_scores[key][2])
	print(key + 'T', exmm_probs[key][3], exmm_scores[key][3])
	print('')


root, ext = args.extxt.split('.')
filename = root + '_mm' + '.tsv'
print(filename)

with open(filename, 'w', newline='') as tsvfile:
	writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
	for key in exmm_scores:
		writer.writerow([key + 'A', exmm_probs[key][0], exmm_scores[key][0]])
		writer.writerow([key + 'C', exmm_probs[key][1], exmm_scores[key][1]])
		writer.writerow([key + 'G', exmm_probs[key][2], exmm_scores[key][2]])
		writer.writerow([key + 'T', exmm_probs[key][3], exmm_scores[key][3]])
tsvfile.close()


#donors = ml.read_txt_seqs(args.dntxt)

#don_pwm, don_ppm = ml.pwm_score(donors)


