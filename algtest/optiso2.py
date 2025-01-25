import argparse
import subprocess
import random
import sys
import json

# ln -s isoform.py into working directory
import isoform
from isoform import Locus

parser = argparse.ArgumentParser()
parser.add_argument('config', type=str, metavar='<json>',
    help='configuration file with all genes')
#parser.add_argument('apc_gen_gene', type=str, metavar='<file>',
#    help='apc generated gff file')
parser.add_argument('--program', required=False, type=str, default='geniso2',
	metavar='<exec>', help='path to program [%(default)s]')

args = parser.parse_args()

# add geniso2 to path
# only need to add this to .bashrc
# export PATH=$PATH:/home/carlos/Code/isoforms

# run this to get a test gff
# ch.13301 is fast to test
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
sys.exit()
'''

# store all isoforms in memory
# don't want to read and store a ton of gffs
# will be much smaller
# fortnite the genetic algorithm
# check if true isoform is in population to start with
# this might ignore important lower prob isoforms
with open(args.config, 'r') as file:
    data = json.load(file)

minin = data['params']['minin']
minex = data['params']['minex']
flank = data['params']['flank']
limit = data['params']['limit']

dpwm = isoform.read_pwm(data['models']['dpwm'])
apwm = isoform.read_pwm(data['models']['apwm'])
elen = isoform.read_len(data['models']['elen'])
ilen = isoform.read_len(data['models']['ilen'])
emm = isoform.read_markov(data['models']['emm'])
imm = isoform.read_markov(data['models']['imm'])

models = (dpwm, apwm, emm, imm, elen, ilen)
models = (None, None, None, None, None, None)
weights = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
weights = (None, None, None, None, None, None)

# converted with log2(p/0.25)
# 0.25 doesn't matter, just need to convert to score
icost = isoform.prob2score(0.001)
icost = 0

for gene in data['genes']:
    fasta = data['apc_dir'] + data['genes'][gene][0]
    gff3 = data['apc_dir'] + data['genes'][gene][1]
    print(gff3)
    name, seq = next(isoform.read_fasta(fasta))
    locus = Locus(name, seq, minin, minex, flank, models, weights, 
                  icost, limit=100)
    isos = locus.isoforms
    for tx in isos:
        print(tx)
        s = 0 
        s += isoform.score_apwm(apwm, tx) * wacc
        s += isoform.score_apwm(apwm, tx) * wacc
        s += isoform.score_apwm(apwm, tx) * wacc
        s += isoform.score_apwm(apwm, tx) * wacc
        s += isoform.score_apwm(apwm, tx) * wacc
        s += isoform.score_apwm(apwm, tx) * wacc



# schema
'''
need a starting population, default 100 like in optiso
each individual is a set of weights
100 sets of weights, each set scores the isoform space
the genes are the weights, and the fitness is the score of the isoform
we are tyring to decrease the manhattan distance
default 50% of population dies after each generation, like in optiso 
'''



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



     

