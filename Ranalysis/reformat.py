import argparse
import os

parser = argparse.ArgumentParser(description="Reformat isorandom output \
    into dataframe format for R")

parser.add_argument('isoran_dir', type=str, metavar='<directory>', \
    help='directory with isorandom output')

args = parser.parse_args()

for fname in os.listdir(args.isoran_dir):
    if fname.endswith('.isorandom'): continue
    splices = fname.split('.')[1]
    with open(f"{args.isoran_dir}{fname}") as fp:
        sims = 0
        time = 0
        for line in fp.readlines():
            if line.startswith('#'): 
                time = line.rstrip().split(' ')[3]
            else:
                sims += 1
        print(splices, sims, time)

'''
rnd.4.950.txt
max splice 4
length 950
seq len, donors, acceptors, isoforms
estimated time
there are 1000 simulations in each file
How to format for use in R?
length, splice, dons, accs, isoforms
length, splice, simulations, elapsed time

isorandom 500 10 --max_splice 2
'''
