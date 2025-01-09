import argparse
import isomod as im
import bt_op as bp
import time

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
                    help='single APC gene to test')
parser.add_argument('apc_dir', type=str, metavar='<directory>',
                    help='directory with APC dataset')

args = parser.parse_args()

# test standard apc algorithm
# ch.13301.fa

seqid, seq = im.read_fasta(args.fasta)

# compare these two
# backtracking algorithm may have issue with first exon position
#seq = 'CCCCCCCGTCCCAGCCCCCGTCCCGTGTCCCAGCCCCCCCAGCCCCCCC'
#seq = 'CCCCCCGTCCCAGCCCCGTGTCCCAGCCCAGCCCCCCC'

seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'


maxs = 3
minin = 3
minex = 3
flank = 3

'''
maxs = 3
minin = 25
minex = 35
flank = 100
'''

sdons, saccs = im.get_gtag(seq, flank, minex)

start_time1 = time.time()
sisos, strials = im.apc(sdons, saccs, maxs, minin, minex, flank, seq)
end_time1 = time.time()

elapsed_time1 = end_time1 - start_time1

print("##### standard #####")
print(elapsed_time1)
print(sdons, saccs)
print(f"isos: {len(sisos)} trials: {strials}")

print("##### backtrack #####")
bdons, baccs = bp.gtag_sites(seq, flank, minex)

start_time2 = time.time()
bisos, info = bp.all_possible(seq, minin, minex, maxs, flank)
end_time2 = time.time()

elapsed_time2 = end_time2 - start_time2
print(elapsed_time2)

print(bdons, baccs)
print(f"isos: {len(bisos)} trials: {info['trails']}")

smRNAs = []
for iso in sisos:
    mRNA = (iso['exons'], iso['introns'])
    smRNAs.append(mRNA)

bmRNAs = []
for iso in bisos:
    mRNA = (iso['exons'], iso['introns'])
    bmRNAs.append(mRNA)
'''
print('#####')
for mRNA in sorted(smRNAs):
    print(mRNA)
print('#####')
for mRNA in sorted(bmRNAs):
    print(mRNA)
'''

# notes
# for smaller genes, backtrack runs more trials (sometimes)
# for larger genes, standard runs more trials
# it is much faster even if there are more trials (re check this)
# large discrepency between number of isoforms with larger genes
# not sure if trials is reported in the same way

# coding a recursive backtracking algorithm with help from youtube

print('##########')

set1 = [6, 17, 22]
set2 = [13, 28]

# test with 100 items freezes my vm
test = [x for x in range(100)]
nums = set1

n = len(nums)
res, sol = [], []

# i for index
def backtrack(i):
    
    # base case, index is at end of list (length n of list)
    if i == n:
        # append copy of sol not reference
        res.append(sol[:])
        return
    
    # left path, don't pick nums[i]
    # move on to next index
    backtrack(i+1)

    # right path, pick nums[i]
    sol.append(nums[i])
    # move on to next leaf
    backtrack(i+1)
    # undo changes, recursively backtrack
    # removes item at last position
    sol.pop()


# start at index 0
#backtrack(0)
print(sorted(res))

r, s = [2], [3]

r.append(s[:])

import itertools

combos = []
for i in range(len(nums)+1):
    for combo in itertools.combinations(nums, i):
        combos.append(combo)

print(sorted(combos))
# backtracking is another way to do all combinations

# itertools will not limit recursion depth, will crash your computer instead
# does this matter in practical applications?

# what is the most amount of donor/acceptor sites in an apc gene?

import os

minex = 35
flank = 100

sdons, saccs = im.get_gtag(seq, flank, minex)

lengths = []
for file in os.listdir(args.apc_dir):
    if file.endswith('gff3'): continue
    else:
        seqid, seq = im.read_fasta(f"{args.apc_dir}{file}")
        dons, accs = im.get_gtag(seq, flank, minex)
        lengths.append(len(dons))
        lengths.append(len(accs))

print(max(lengths))

# at most there are 95 dons/accs sites in one gene
# so don't need to worry about hitting recursion depth with backtracking
# my vm with 8 cores freezes when i track backtracking with a list of 100
# itertools runs just fine
# not sure why the dif, need to limit cores used
# not sure how many gtag sites Gong tested

