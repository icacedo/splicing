import argparse
import isoform

if __name__ == '__main__':

	## Command Line Interface ##

	parser = argparse.ArgumentParser(
		description='Alternative Splice Generator')
	parser.add_argument('fasta', type=str, metavar='<file>', help='fasta file')
	parser.add_argument('--intron', required=False, type=int, default=35,
		metavar='<int>', help='minimum length of intron [%(default)i]')
	parser.add_argument('--exon', required=False, type=int, default=25,
		metavar='<int>', help='minimum length exon [%(default)i]')
	parser.add_argument('--splice', required=False, type=int, default=3,
		metavar='<int>', help='maximum number of introns [%(default)i]')
	parser.add_argument('--flank', required=False, type=int, default=99,
		metavar='<int>', help='distance to ignore on each side [%(default)i]')
	parser.add_argument('--dpwm', required=False, type=str, metavar='<file>',
		help='position weight matrix for donor site [%(default)s]')
	parser.add_argument('--apwm', required=False, type=str, metavar='<file>',
		help='position weight matrix for acceptor site [%(default)s]')
	parser.add_argument('--emm', required=False, type=str, metavar='<file>',
		help='markov model for exon sequence [%(default)s]')
	parser.add_argument('--imm', required=False, type=str, metavar='<file>',
		help='markov model for intron sequence [%(default)s]')
	parser.add_argument('--elen', required=False, type=str, metavar='<file>',
		help='length model for exons [%(default)s]')
	parser.add_argument('--ilen', required=False, type=str, metavar='<file>',
		help='length model for introns [%(default)s]')
	parser.add_argument('--gff', required=False, type=str, metavar='<file>',
		help='use GFF for source of splice sites [%(default)s]')
	parser.add_argument('--full', action='store_true',
		help='see all splice forms')
	parser.add_argument('--limit', required=False, type=int, metavar='<int>',
		help='limit full report')
	arg = parser.parse_args()

	for id, seq in isoform.read_fasta(arg.fasta):
		txs, info = isoform.all_possible(seq, arg.intron, arg.exon,
			arg.splice, arg.flank, gff=arg.gff)
		print('seq:', id)
		print('len:', len(seq))
		print('donors:', info['donors'])
		print('acceptors:', info['acceptors'])
		print('trials:', info['trials'])
		print('isoforms:', len(txs))

		if arg.full:
			deets = arg.limit if arg.limit else len(txs)
			print('details:', deets)

			# sort by score first

			for i in range(deets):
				tx = txs[i]
				for intron in tx['introns']:
					print(intron[0], '..', intron[1], ' ', sep='', end = '')
				print()

				# output might rather be GFF or something

