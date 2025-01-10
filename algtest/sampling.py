import random
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('apc_dir', type=str, help='directory with APC genes')

args = parser.parse_args()

'''
introns = []
count = 0
for i in range(10, 40, 10):
    int_list = [count for x in range(i)]
    for num in int_list:
        introns.append(num)
    count += 1
   
random.shuffle(introns)

sub_ints = []
for i in range(15):
    sub_ints.append(random.choice(introns))
'''

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

n_samples = 10
increment = 1
start = 1

all_data = {}
for gene in genes:
    data = {}
    for i in range(start, n_samples+start, increment):
        pop = []
        probs = []
        for intron in genes[gene]:
            pop.append(intron)
            probs.append(genes[gene][intron])
       
        # i don't understand how cumulitive weights work 
        sample = random.choices(pop, weights=probs, k=i)

        sampled_introns = {}
        for intron in sample:
            if intron not in sampled_introns:
                sampled_introns[intron] = 1
            else:
                sampled_introns[intron] += 1
        data[i] = sampled_introns
    all_data[gene] = data
    break

# write to csv file for import to R
print(all_data)
for i in all_data:
    for j in all_data[i]:
        print(j)




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


