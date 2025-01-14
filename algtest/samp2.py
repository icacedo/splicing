import sys

# pick gene with most amount of introns
# 15294
apc_gene = sys.argv[1]


pop = []
probs = []
with open(sys.argv[1]) as fp:
    for line in fp.readlines():
        line = line.rstrip()
        line = line.split('\t')
        if line[1] == 'RNASeq_splice' and line[2] == 'intron':
            intron = (line[3], line[4])
            score = line[5]
            print(intron, score)