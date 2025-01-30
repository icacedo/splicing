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
parser.add_argument('--apc_prog', required=False, type=str, default='geniso2',
	metavar='<exec>', help='path to apc program [%(default)s]')
parser.add_argument('--cmp_prog', required=False, type=str, default='cmpiso',
    metavar='<exec>', help='path to isoform comparison program \
        [%(default)s]')

args = parser.parse_args()

# add geniso2 and cmpiso to path
# only need to add this to .bashrc
# export PATH=$PATH:/home/carlos/Code/isoforms

# run this to get a test gff
# ch.13301 is fast to test
'''
subprocess.run([
    'geniso2', '../../datacore2024/project_splicing/smallgenes/ch.2_1.fa',
    '--dpwm', '../../isoforms/models/don.pwm',
    '--apwm', '../../isoforms/models/acc.pwm',
    '--emm', '../../isoforms/models/exon.mm', 
    '--imm', '../../isoforms/models/intron.mm', 
    '--elen', '../../isoforms/models/exon.len',
    '--ilen', '../../isoforms/models/intron.len',
    '--limit', '2'
], stdout=open('ch.2_1.geniso2.gff', 'w'))
sys.exit()
'''

# store all isoforms in memory
# don't want to read and store a ton of gffs
# will be much smaller
# fortnite the genetic algorithm
# check if true isoform is in population to start with
# this might ignore important lower prob isoforms

def random_start():

    single = {
            'genotype': {
                'wdpwm': random.random(),
                'wapwm': random.random(),
                'wemm': random.random(),
                'wimm': random.random(),
                'welen': random.random(),
                'wilen': random.random(),
                'icost': random.random()
            },
            'fitness' : None
        }
    
    return single

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

single = random_start()

# make sure that intron probabilities for cmpiso are used in the
# same way that uses input gffs
command = (
    f'geniso2 ../../datacore2024/project_splicing/smallgenes/ch.2_1.fa '
    f'--dpwm ../../isoforms/models/don.pwm '
    f'--apwm ../../isoforms/models/acc.pwm '
    f'--emm ../../isoforms/models/exon.mm '
    f'--imm ../../isoforms/models/intron.mm '
    f'--elen ../../isoforms/models/exon.len '
    f'--ilen ../../isoforms/models/intron.len '
    f'--wdpwm {single['genotype']['wdpwm']} '
    f'--wapwm {single['genotype']['wapwm']} '
    f'--wemm {single['genotype']['wemm']} '
    f'--wimm {single['genotype']['wimm']} '
    f'--welen {single['genotype']['welen']} '
    f'--wilen {single['genotype']['wilen']} '
    f'--icost {single['genotype']['icost']} '
    f'--limit 3'
)

res = subprocess.run([command], shell=True, text=True, 
                     stdout=open('ch.2_1.geniso2.gff', 'w'))

sys.exit()




# converted with log2(p/0.25)
# 0.25 doesn't matter, just need to convert to score
icost = isoform.prob2score(0.001)
icost = 0

for gene in data['genes']:
    fasta = data['apc_dir'] + data['genes'][gene][0]
    gff3 = data['apc_dir'] + data['genes'][gene][1]
    name, seq = next(isoform.read_fasta(fasta))
    models = (None, None, None, None, None, None)
    weights = (1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
    locus = Locus(name, seq, minin, minex, flank, models, weights, 
                  icost, limit=2)
    isos = locus.isoforms
    #single = random_start()
    weight = []
    total = 0
    for tx in isos:
        s = 0 
        s += isoform.score_apwm(apwm, tx) * single['genotype']['wdpwm']
        s += isoform.score_apwm(apwm, tx) * single['genotype']['wapwm']
        s += isoform.score_apwm(apwm, tx) * single['genotype']['wemm']
        s += isoform.score_apwm(apwm, tx) * single['genotype']['wimm']
        s += isoform.score_apwm(apwm, tx) * single['genotype']['welen']
        s += isoform.score_apwm(apwm, tx) * single['genotype']['wilen']
        s += len(tx['introns']) * single['genotype']['icost']
        tx['score'] = s
        w = 2 ** s
        weight.append(w)
        total += w
        tx['prob'] = w
    for tx in isos:
        tx['prob'] = tx['prob'] / total
    print('###')
    print(isos)
    print('###')
    for x in isos:
        for intron in x['introns']:
            print(intron, x['prob'])
        


# schema
'''
starting population of 100 singles
for each gene, there X number of isoforms
fitness is the manhattan distance between apc and wormbase introns
need to write cmpiso into this script so i don't need to read in two gffs
'''
# ??????????????????????????????????????????????????????????????
# for apc generated isoforms with multiple introns
# the intron score for each intron is the same as the mRNA score

print('##########')
def mdist(apc_gff, g2):

    i1 = isoform.get_introns(apc_gff)
    i2 = isoform.get_introns(g2)

    print(i1, '@')
    print(i2, '@')

mdist(gff3, 'ch.2_1.geniso2.gff')

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
    
    return solo


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



