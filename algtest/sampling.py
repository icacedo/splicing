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

#
                    


