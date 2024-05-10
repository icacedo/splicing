import argparse
import os
import json

parser = argparse.ArgumentParser()
parser.add_argument('wb_dir', type=str, metavar='<directory>',
	help='path to directory with wormbase fasta and gff files')
parser.add_argument('model_dir', type=str, metavar='<directory>',
	help='path to directory with apc model files')
parser.add_argument('--gff', required=False, type=str, 
	help='use wormbase gff donor/acceptor sites')

parser.add_argument('--min_intron', required=False, type=int, default=35,
	metavar='<int>', help='minimum length of intron [%(default)i]')
parser.add_argument('--min_exon', required=False, type=int, default=25,
	metavar='<int>', help='minimum length exon [%(default)i]')
parser.add_argument('--max_splice', required=False, type=int, default=3,
	metavar='<int>', help='maximum number of introns [%(default)i]')
parser.add_argument('--flank', required=False, type=int, default=100,
	metavar='<int>', help='genomic flank on each side [%(default)i]')
parser.add_argument('--limit', required=False, type=int, default=100,
	metavar='<int>', help='limit number of transcripts [%(default)i]')

arg = parser.parse_args()

config = {}

if arg.gff: config['gff_introns'] = 'true'
else: config['gff_introns'] = 'false'

apwm = None
dpwm = None
emm = None
imm = None
elen = None
ilen = None
for file in os.listdir(arg.model_dir):
	if file == 'acc.pwm': apwm = arg.model_dir + file
	if file == 'don.pwm': dpwm = arg.model_dir + file
	if file == 'exon.mm': emm = arg.model_dir + file
	if file == 'intron.mm': imm = arg.model_dir + file
	if file == 'exon.len': elen = arg.model_dir + file
	if file == 'intron.len': ilen = arg.model_dir + file

config['cli'] = {
	'--min_exon': arg.min_exon,
	'--min_intron': arg.min_intron,
	'--max_splice': arg.max_splice,
	'--flank': arg.flank,
	'--limit': arg.limit,
	'--apwm': apwm,
	'--dpwm': dpwm,
	'--emm': emm,
	'--imm': imm,
	'--elen': elen,
	'--ilen': ilen
}


fpaths = {}
for file in os.listdir(arg.wb_dir):
	fsp = file.split('.')
	name = f'{fsp[0]}.{fsp[1]}' 
	if name not in fpaths:
		fpaths[name] = []
		if file.endswith('.fa'): fpaths[name].append(arg.wb_dir+file)
		if file.endswith('.gff3'): fpaths[name].append(arg.wb_dir+file)
	else:
		if file.endswith('.fa'): fpaths[name].append(arg.wb_dir+file)
		if file.endswith('.gff3'): fpaths[name].append(arg.wb_dir+file)

fpaths = {name:sorted(fpaths[name]) for name in fpaths.keys()}

def add_gene(config, name, paths):

	gconfig = config
	ginfo = {
		'name': name,
		'fasta': paths[0],
		'gff': paths[1]
	}
	gconfig['data'] = [ginfo]

	return gconfig

c1 = 0
c2 = 1
for gene in fpaths:
	gconfig = add_gene(config, gene, fpaths[gene])
	iid = gene.split('.')[1]
	if not os.path.exists(f'outfigs{c2}/'): os.makedirs(f'outfigs{c2}/')
	with open(f'outfigs{c2}/{iid}.config.json', 'w') as jfile:
		json.dump(gconfig, jfile, indent=4)
	c1 += 1
	if c1%70 == 0: c2 += 1




