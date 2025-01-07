import argparse
import isomod as im
import bt_op as bp
import time

parser = argparse.ArgumentParser()
parser.add_argument('fasta', type=str, metavar='<file>',
                    help='single APC gene to test')

args = parser.parse_args()

# test standard apc algorithm
# ch.13301.fa

seqid, seq = im.read_fasta(args.fasta)

# compare these two
# backtracking algorithm may have issue with first exon position
seq = 'CCCCCCCGTCCCAGCCCCCGTCCCGTGTCCCAGCCCCCCCAGCCCCCCC'
seq = 'CCCCCCGTCCCAGCCCCGTGTCCCAGCCCAGCCCCCCC'

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


don = []
acc = []
isoforms = []

def backtrack(i):

    if i == maxs * 2: return

    if i % 2 == 0:

        for ds in bdons:
            if don and ds <= acc[-1]: continue
            else:
                if ds - 1 - flank < minex:
                    continue
            if acc and ds < acc[-1] + minex + 2:
                continue

            don.append(ds)
            backtrack(i + 1)
            don.pop()
        
    else:

        for ac in acc:

            if ac <= don[-1]: continue
            if ac < don[-1] + minin - 1:
                continue

            acc.append(ac)

            if len(seq) - flank - ac + 1 >= minex:
                tx = bp.build_mRNA(seq, flank, len(seq) - flank - 1, don, acc)
                isoforms.append(tx)
                backtrack(i + 1)
                acc.pop()
            else:
                acc.pop()
                continue

backtrack(0)

print(isoforms, info)
            




