import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('jb_gff')

args = parser.parse_args()

lines = ['wow', 'wowow']

def writer(lines, outfile, feature):
	
	fpath = f'{outfile}.{feature}.gff3'
	if os.path.exists(fpath):
		os.remove(fpath)
	f = open(fpath, 'a')
	for line in lines:
		f.write(f'{line}\n')
	f.close()

isos = []
ints = []
exos = []
dons = []
accs = []
with open(args.jb_gff, 'r') as fp:
	count = 0 
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('#'):
			count += 1
			continue
		if count == 1:
			isos.append(line)
		if count == 2:
			ints.append(line)
		if count == 3:
			exos.append(line)
		if count == 4:
			dons.append(line)
		if count == 5:	
			accs.append(line)

bname = os.path.splitext(os.path.basename(args.jb_gff))[0]
writer(isos, bname, 'isos')
writer(ints, bname, 'ints')
writer(exos, bname, 'exos')
writer(dons, bname, 'dons')
writer(accs, bname, 'accs')
			
			
			
