import random
import argparse
import os
import csv

parser = argparse.ArgumentParser()
parser.add_argument('apc_dir', type=str, help='directory with APC genes')
parser.add_argument('--n_samples', type=int, 
                    help='choose how many times a population is sampled')
parser.add_argument('--increment', type=int, help='increment sampling number')
parser.add_argument('--start', type=int, help='sample starting size')

args = parser.parse_args()

genes = {}
for fname in os.listdir(args.apc_dir):
    if fname.endswith('gff3'):
        gid = fname.split('.')[1]
        with open(f'{args.apc_dir}{fname}') as fp:
            intron_counts = {}
            for line in fp.readlines():
                line = line.rstrip()
                line = line.split('\t')
                if line[1] == 'RNASeq_splice' and line[2] == 'intron':
                    intron = (line[3], line[4])
                    reads = line[5]
                    intron_counts[intron] = reads
            genes[gid] = intron_counts

for gene in genes:
    total = 0
    for intron in genes[gene]:
        total += float(genes[gene][intron])
    for intron in genes[gene]:
        genes[gene][intron] = float(genes[gene][intron])/total
        #print(gene, intron, genes[gene][intron])

# each intron for each gene now has a probability

n_samples = args.n_samples
increment = args.increment
start = args.start

gene_samples = {}
for g in genes:
    pop = []
    probs = []
    for i in genes[g]:
        pop.append(i)
        probs.append(genes[g][i])
    sampled_introns = []
    for j in range(start, n_samples+increment, increment):
        sample = random.choices(pop, weights=probs, k=j)
        print(f'{g} sampling done')
        sampled_introns.append(sample)
    gene_samples[g] = sampled_introns

data = {}
for g in gene_samples:
    run_info = {}
    for int_group in gene_samples[g]:
        sums = {}
        for intron in int_group:
            if intron not in sums:
                sums[intron] = 1
            else:
                sums[intron] += 1
        run_info[len(int_group)] = sums
    data[g] = run_info

# write to csv file for import into R
with open('depth_sim.csv', 'w') as csvfile:
    datawriter = csv.writer(csvfile)
    datawriter.writerow(['gene_id', 'sim_num', 'intron', 'icount'])
    for gene in data:
        for info in data[gene]:
            for intron in data[gene][info]:
                datawriter.writerow([gene, info, 
                                    f'b{int(intron[0])}e{int(intron[1])}', 
                                    data[gene][info][intron]])

     








'''
pop = [1, 2, 3, 4]
w = [0.1, 0.3, 0.2, 0.4]

sums = {}
for i in range(1000):
    num = random.choices(pop, weights=w, k=1)[0]
    #nummers = random.choices(pop, weights=w, k=3)
    if num not in sums:
        sums[num] = 1
    else:
        sums[num] += 1

total = sum(sums.values())
probs = {key: num / total for key, num in sums.items()}

# sampled probs are roughly equal 
# to weights
print(dict(sorted(probs.items())))

numbers = random.choices(pop, weights=w, k=1000)

all = {}
for n in numbers:
    if n not in all:
        all[n] = 1
    else:
        all[n] += 1

print(all)
tot = sum(all.values())
print(tot)
ps = {key: n / tot for key, n in all.items()}
# everything looks correct
print(dict(sorted(ps.items())))
'''


