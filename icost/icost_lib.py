import pickle
import modelib as ml
import sys

fasta = sys.argv[1]

def apc_pickler(fasta_path, maxs=3, min_intron=25, 
				min_exon=25, flank=100, gff=None):

	seqid = None
	seq = None
	for seqid, seq in ml.read_fastas(fasta):
		seqid = seqid
		seq = seq

	if gff: 
		dons, accs = ml.read_gff_sites(seq, gff)
	else:
		dons, accs = ml.get_gtag(seq)
		
	apc_isoforms, trials = ml.apc(dons, accs, maxs, minin, minex, flank, seq)

	fname = fasta_path.split('/')[-1]
	ID = fname.split('.')[1]
	outdir = args.outdir+'apc_pickles/'
	os.makedirs(os.path.dirname(outdir), exist_ok=True)

	if gff:
		name = outdir+'ch.'+ID+'.apc_isoforms.bli.pkl'
	else:
		name = oudir+'ch.'+ID+'.apc_isoforms.pkl'

	with open(name, 'rb') as pick:
		pickle.dump(apc_isoforms, pick)

	with open(name, 'rb') as pick:
		pickell = pickle.load(pick)

	assert apc_isoforms == pickell, 'pickled incorrectly'





wow = 4
apc_pickler(fasta, 'wow', wow)
	
