import random
import statistics
import argparse
import os
import csv

# pick gene with most amount of introns
# 15294

parser = argparse.ArgumentParser()
parser.add_argument('apc_dir', type=str, help='directory with APC genes')
parser.add_argument('--n_experiments', type=int,
                    help='choose how many experiments to run')
parser.add_argument('--n_samples', type=int, 
                    help='choose how many times a population is sampled')
parser.add_argument('--increment', type=int, help='increment sampling number')
parser.add_argument('--start', type=int, 
                    help='sample starting size, start should equal increment')

args = parser.parse_args()

path = args.apc_dir
files = [[path, file] for file in os.listdir(path) if file.endswith('gff3')]

all_genes = {}
c = 0
for path, file in files:
    pop = []
    scores = []
    #with open(f'{path}{file}') as fp:
    with open(f'{path}ch.15294.gff3') as fp:
        gid = file.split('.')[1]
        for line in fp.readlines():
            line = line.rstrip()
            line = line.split('\t')
            if line[1] == 'RNASeq_splice' and line[2] == 'intron':
                pop.append((line[3], line[4]))
                scores.append(float(line[5]))

    icount = {}
    for exp in range(args.n_experiments):
        seen = set()
        for i in range(args.start, args.n_samples+args.start, args.increment):
            for intron in random.choices(pop, weights=scores, k=i):
                seen.add(intron)
            n = len(seen)
            if i not in icount: icount[i] = []
            icount[i].append(n)

    means = {}
    for samples, data in icount.items():
        means[samples] = statistics.mean(data)
        #print(samples, statistics.mean(data), flush=True)

    all_genes[gid] = means
    
    c += 1
    if c == 2: break

for i in all_genes.items():
    print(i)

print('#####')








# don't need to calculate the all gene averages in R
'''
col_names = [i for i in range(args.start, args.n_samples+args.start, 
                              args.increment)]

col_names.insert(0, 'gid')

print(col_names)


for g in all_genes:
    print(g, f'{[n_samp[1] for n_samp in all_genes[g].items()]}')

print(f'{[n for n in range(10)]}')

col_names = [i for i in range(10, 100, 10)]
col_names.insert(0, 'gid')
print(col_names)
'''

'''
with open('samp_sim.csv', 'w') as csvfile:
    datawriter = csv.writer(csvfile)
    datawriter.writerow(['gene_id', 'sim_num', 'intron', 'icount'])
    for gene in data:
        for info in data[gene]:
            for intron in data[gene][info]:
                datawriter.writerow([gene, info, 
                                    # each intron label needs to be unique
                                    f'{gene}b{int(intron[0])}e{int(intron[1])}', 
                                    data[gene][info][intron]])
'''                

'''
# test scores vs probs

items = [1, 2, 3, 4]
probs = [0.1, 0.3, 0.2, 0.4]
scores = [10, 30, 20, 40]

n = 100

picks1 = []
picks2 = []
for pick in random.choices(items, weights=probs, k=n):
    picks1.append(pick)
for pick in random.choices(items, weights=scores, k=n):
    picks2.append(pick)

summary1 = {}
for p in picks1:
    if p not in summary1:
        summary1[p] = 1
    else:
        summary1[p] += 1

summary2 = {}
for p in picks2:
    if p not in summary2:
        summary2[p] = 1
    else:
        summary2[p] += 1

s1 = {key: summary1[key] for key in sorted(summary1)}
s2 = {key: summary2[key] for key in sorted(summary2)}

print(s1)
print(s2)
'''