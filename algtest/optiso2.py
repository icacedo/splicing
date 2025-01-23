import argparse
import subprocess
import random

parser = argparse.ArgumentParser()
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



     

