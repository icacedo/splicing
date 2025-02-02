import isoform2

seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'

maxs = 3
minin = 3
minex = 3
flank = 3

dons, accs = isoform2.gtag_sites(seq, flank, minex)

print(dons, accs)

introns = []
for d in dons:
    for a in accs:
        introns.append((d, a))

print(introns)

print('#####')

import copy

# https://medium.com/algorithms-and-leetcode/backtracking-e001561b9f28

def perm(a, n, k, depth, used, curr, ans):

    if depth == k:
        #ans.append(curr[::])
        ans.append(copy.copy(curr))
        return
    
    for i in range(n):
        if not used[i]:
            curr.append(a[i])
            used[i] = True
            #print(curr)
            perm(a, n, k, depth+1, used, curr, ans)

            curr.pop()
            #print('backtrack: ', curr)
            used[i] = False
    return

a = [1, 2, 3]
n = len(a)
ans = [[None]]
used = [False] * len(a)
ans = []

perm(a, n, n, 0, used, [], ans)
print(ans)

n = len(introns)
ans =[[None]]
used = [False] * len(introns)
ans = []

dons = [6, 18, 24]
accs = [12, 30]

introns = []
for d in dons:
    for a in accs:
        if d > a: continue
        intron = (d, a)
        introns.append(intron)

print(introns)


bad_iso = [(24, 30), (18, 30), (6, 12), (6, 30)]

def check_iso(bad_iso):

    for intron in bad_iso:
        if intron == bad_iso[0]:
            print('first', intron)

        if intron == bad_iso[-1]:
            print('last', intron)

check_iso(bad_iso)












