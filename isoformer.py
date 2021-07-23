import argparse
import itertools
import korflib

loopcount = 0

def isoforms(don, acc, imin, emin):
	global loopcount
	for d in range(0, len(don)):
		for dsites in itertools.combinations(don, d):
			for a in range(0, len(acc)):
				for asites in itertools.combinations(acc, a):
					loopcount += 1
		
					# same number of sites both sides required
					if len(dsites) != len(asites): continue
			
					# make transcript based on intron positions
					introns = []
					for i in range(len(dsites)):
						if dsites[i] < asites[i]:
							introns.append((dsites[i], asites[i]))
					if len(introns) == 0: continue
			
					# check intron lengths
					introns_ok = True
					for site in introns:
						l = site[1] - site[0] + 1
						if l < imin:
							introns_ok = False # too short
							break
					if not introns_ok: continue
					
					# check exon lengths
					exons_ok = True
					for i in range(1, len(introns)):
						beg = introns[i-1][1] +1
						end = introns[i][0] -1
						l = end - beg + 1
						if l < emin:
							exons_ok = False
							break
					if not exons_ok: continue
					
					# done with this isoform (which might not be unique)
					yield introns


if __name__ == '__main__':

	## Command Line Interface ##

	parser = argparse.ArgumentParser(
		description='Theoretical isoform enumerator')
	parser.add_argument('fasta', type=str, metavar='<file>', help='fasta file')                
	parser.add_argument('--intron', required=False, type=int, default=35,
		metavar='<int>', help='minimum length of intron [%(default)i]')
	parser.add_argument('--exon', required=False, type=int, default=25,
		metavar='<int>', help='minimum length exon [%(default)i]')
	parser.add_argument('--full', action='store_true',
		help='see all splice forms')
	arg = parser.parse_args()

	for id, seq in korflib.read_fasta(arg.fasta):
		# find donor sites
		don = []
		for i in range(arg.exon, len(seq) -arg.exon):
			if seq[i:i+2] == 'GT': don.append(i)
	
		# find acceptor sites
		acc = []
		for i in range(arg.exon, len(seq) -arg.exon):
			if seq[i:i+2] == 'AG': acc.append(i+1)
	
		# find non-redundant isoforms
		isoform = {}
		isocount = 0
		for tx in isoforms(don, acc, arg.intron, arg.exon):
			txs = f'{tx}'
			isocount += 1
			if txs not in isoform:
				isoform[txs] = tx
	
		# output	
		print(id, len(isoform), isocount, loopcount)
		if arg.full:
			for iso in isoform:
				print(iso)
	
