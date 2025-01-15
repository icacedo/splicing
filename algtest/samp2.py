import sys
import random
import statistics

# pick gene with most amount of introns
# 15294
apc_gene = sys.argv[1]


pop = []
scores = []
with open(sys.argv[1]) as fp:
    for line in fp.readlines():
        line = line.rstrip()
        line = line.split('\t')
        if line[1] == 'RNASeq_splice' and line[2] == 'intron':
            pop.append((line[3], line[4]))
            scores.append(float(line[5]))

# 

icount = {}
for exp in range(5):
    seen = set()
    for i in range(10, 1000, 10):
        for intron in random.choices(pop, weights=scores, k=i):
            seen.add(intron)
        n = len(seen)
        if i not in icount: icount[i] = []
        icount[i].append(n)

for samples, data in icount.items():
    print(samples, statistics.mean(data), flush=True)