import argparse
import json
import os
import isosort_lib as isl

parser = argparse.ArgumentParser()
parser.add_argument('json', type=str, metavar='<file>',
	help='single gene json file')
parser.add_argument('fasta', type=str, metavar='<file>',
	help='single gene fasta file')

args = parser.parse_args()
	
bname = os.path.basename(args.json)
gID = bname.split('.')[0]

seq = isl.get_seq(args.fasta)
print(len(seq))

with open(args.json, 'r') as jfile:
	info = json.load(jfile)
	print(info[f'ch.{gID}-1'])
	wbgene = info[f'ch.{gID}-wb']['Parent=Gene']
	
# i think viewbox is based on gene size
with open('test.html', 'w') as hfile:
	hfile.write(f'<html><title>{gID}</title><body>\n')
	hfile.write(f'<h1>ch.{gID} {wbgene}</h1>\n')
	hfile.write('<svg viewBox="0 0 902 200" xmlns="http://www.w3.org/2000/svg">\n')	
	hfile.write('<rect fill="#ddd" x="0" y="0" width="100%" height="100%"/>\n')
	hfile.write('<rect x="0" y="40" width="100" height="12" fill="#949494"></rect>\n')
	hfile.write('<rect x="100" y="40" width="113" height="12" fill="#3368FF"></rect>\n')
	hfile.write('<rect x="200" y="40" width="220" height="12" fill="#FFD433"></rect>\n')
	hfile.write('<font>')
	hfile.write('<font face\n')
	hfile.write('font-family="Super Sans"
	hfile.write('
	hfile.write('</font>')
	hfile.write('<text x="20" y="35" class="sans">ACGT</text>\n')
	hfile.write('</svg>\n')
	hfile.write('</body></html>\n')
	















		


