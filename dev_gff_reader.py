import sys
import modelib
import isoform

# get gt/ag sites from gff file
# use ch.9940

def read_fastas(fastafile):

	with open(fastafile) as fp:

		seqid = ''
		seq = ''	

		for line in fp.readlines():
			line = line.rstrip()
			if line.startswith('>'): seqid += line
			else: seq += line
		yield (seqid, seq)			
	fp.close()		

def gff_sites(seq, gff, gtag=True):

	dons = []
	accs = []

	with open(gff) as fp:
		while True:
			line = fp.readline()
			line = line.rstrip()
			if not line: break
			fields = line.split('\t')
			if fields[2] == 'intron':
				beg = int(fields[3]) - 1
				end = int(fields[4]) - 1
				if gtag:
					if seq[beg:beg+2] == 'GT': dons.append(beg)
					if seq[end-1:end+1] == 'AG': accs.append(end)
	fp.close()	
	return sorted(set(dons)), sorted(set(accs))

for i in read_fastas(sys.argv[2]):
	dons, accs = gff_sites(i[1], sys.argv[1])
	print(dons, accs)
	dons1, accs1 = modelib.get_gtag(i[1])
	print(dons1, accs1)	
	dons2, accs2 = isoform.gff_sites(i[1], sys.argv[1])
	print(dons2, accs2)	



