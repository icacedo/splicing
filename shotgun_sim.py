'''
How much sequencing do you have to do to ensure that you get enough 
coverage to examine someone's genotype? Make a simulation of shotgun
sequencing where the variables include:
-size of the genome
-length of the sequencing read
-intended coverage of the genome
Your simulation should "align" reads onto the genome and then examine 
depth of coverage after alignment. Make a table showing intended depth
of coverage vs. actual depth of coverage for various intended levels of 
coverage. For example, you might intend 4x, but some areas will be 1x, 
or even 0x.
-How does read size or genome size affect the simulation?
-How much sequencing does your company propose to perform?
-Is there another way you could have performed this experiment?
'''

import argparse
import random

parser = argparse.ArgumentParser(
	description='sequencing coverage simulator')
parser.add_argument('--genome_size', type=int, metavar='<integer>', 
	default=100, help='size of the genome [%(default)i]')
parser.add_argument('--read_len', type=int, metavar='<integer>',	
	default=10, help='length of the sequencing read [%(default)i]')
parser.add_argument('--coverage', type=int, metavar='<integer>',	
	default=4, help='intended coverage of the genome [%(default)i]')

args = parser.parse_args()

gs = args.genome_size
rl = args.read_len
ce = args.coverage

random.seed()

frag1 = (1, gs)
		
def get_frags(frag1, rl):

	frags = []
	def cutter(frag1, rl):
		
		f_len = frag1[1] - frag1[0] + 1
		if f_len >= rl * 2:
			cut = random.randint(frag1[0], frag1[1])
			frag2 = (frag1[0], cut)
			frag3 = (cut+1, frag1[1])
			cutter(frag2, rl)
			cutter(frag3, rl)
		else:
			frags.append(frag1)
		
	cutter(frag1, rl)
	
	return frags

multi_gens = []
f_lens = []
for i in range(ce):
	frags = get_frags(frag1, rl)
	multi_gens.append(frags)
	for f in frags:
		f_len = f[1] - f[0] + 1
	f_lens.append(f_len)
	
print(sum(f_lens)/len(f_lens))
	
#for m in multi_gens:
#	print(m)





