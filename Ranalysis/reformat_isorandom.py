import argparse
import os
import csv

parser = argparse.ArgumentParser(description="Reformat isorandom output \
    into dataframe format for R")

parser.add_argument('isoran_dir', type=str, metavar='<directory>', \
    help='directory with isorandom output')

args = parser.parse_args()

by_isos = []
by_times = []
for fname in os.listdir(args.isoran_dir):
    if fname.endswith('.isorandom'): continue
    splices = fname.split('.')[1]
    length = fname.split('.')[2]
    with open(f"{args.isoran_dir}{fname}") as fp:
        sims = 0
        time = 0
        for line in fp.readlines():
            line = line.rstrip()
            if line.startswith('#'): 
                time = line.split(' ')[3]
            else:
                sims += 1
                values = [int(value) for value in line.split('\t')]
                values.append(int(splices))
                by_isos.append(values)
        by_times.append([int(length), int(splices), int(sims), float(time)])

header_t = ['length', 'introns', 'simulations', 'seconds']
with open('isoran_times.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(header_t)
    writer.writerows(by_times)

header_i = ['length', 'donors', 'acceptors', 'isoforms', 'introns']
with open('isoran_isos.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(header_i)
    writer.writerows(by_isos)
