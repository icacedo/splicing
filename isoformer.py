import argparse
import korflib
import isoform



if __name__ == '__main__':

	## Command Line Interface ##

	parser = argparse.ArgumentParser(
		description='Theoretical isoform enumerator')
	parser.add_argument('fasta', type=str, metavar='<file>', help='fasta file')                
	parser.add_argument('--intron', required=False, type=int, default=35,
		metavar='<int>', help='minimum length of intron [%(default)i]')
	parser.add_argument('--exon', required=False, type=int, default=25,
		metavar='<int>', help='minimum length exon [%(default)i]')
	parser.add_argument('--max', required=False, type=int, default=3,
		metavar='<int>', help='maximum number of introns [%(default)i]')
	parser.add_argument('--full', action='store_true',
		help='see all splice forms')
	arg = parser.parse_args()

	for id, seq in korflib.read_fasta(arg.fasta):
		txs, info = isoform.all_possible(seq, arg.intron, arg.exon, arg.max)
		print(id, len(txs), info)
		if arg.full:
			for tx in txs:
				for intron in tx:
					print(intron['beg'], '..', intron['end'], ' ', sep='', end = '')
				print()
		
