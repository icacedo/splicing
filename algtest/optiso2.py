import argparse
import subprocess
import random
import json

parser = argparse.ArgumentParser()
#parser.add_argument('config', required=True, type=str, metavar='<json>',
#    help='configuration file with all genes')
parser.add_argument('apc_gen_gene', type=str, metavar='<file>',
    help='apc generated gff file')
parser.add_argument('--program', required=False, type=str, default='geniso2',
	metavar='<exec>', help='path to program [%(default)s]')

arg = parser.parse_args()

# add geniso2 to path
# only need to add this to .bashrc
# export PATH=$PATH:/home/carlos/Code/isoforms

# schema
'''
need a starting population, default 100 like in optiso
what is an individual?
the genes are the weights, and the fitness is the score of the isoform
we are tyring to decrease the manhattan distance
default 50% of population dies after each generation, like in optiso 
'''

# ch.13301 is fast to test
seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'

maxs = 3
minin = 3
minex = 3
flank = 3

# generate isoforms

subprocess.run(['ls'])

popn = 10

def random_start():

    solo = {
            'genotype': {
                '--wdpwm': random.random(),
                '--wapwm': random.random(),
                '--wemm': random.random(),
                '--wimm': random.random(),
                '--welen': random.random(),
                '--wilen': random.random(),
                '--icost': random.random()
            },
            'fitness' : None
        }

# since all isoforms will be generated before hand
# geniso2 does not need to be called here
# input a list with all geniso2 output

# run this to get a test gff
'''
subprocess.run([
    'geniso2', '../../isoforms/apc/ch.13301.fa',
    '--dpwm', '../../isoforms/models/don.pwm',
    '--apwm', '../../isoforms/models/acc.pwm',
    '--emm', '../../isoforms/models/exon.mm', 
    '--imm', '../../isoforms/models/intron.mm', 
    '--elen', '../../isoforms/models/exon.len',
    '--ilen', '../../isoforms/models/intron.len',
], stdout=open('ch.13301.geniso2.gff', 'w'))
'''

'''
class Evo:

    def __init__(self, solo):
        self.solo = {
            'genotype': {
                '--wdpwm': random.random(),
                '--wapwm': random.random(),
                '--wemm': random.random(),
                '--wimm': random.random(),
                '--welen': random.random(),
                '--wilen': random.random(),
                '--icost': random.random()
            },
            'fitness' : None
        }


guy = Evo()

print(guy.solo)
'''
# pseudo
'''
pop1 = class_stuff(init)
pop2 = class_stuff(pop1)
pop3 = class_stuff(pop2)

evo.best_i = best set of weights for the individual
evo.
'''



     

