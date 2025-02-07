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




############
'''
Don't need models, weights to generate isoforms
self.exons/introns = [coordinates]
transcript class just computes score
Locus class contains isoforms
self.worst whatever is currently the worst isoform is stored, skips over bad stuff
this just tells code to end if you are over the isoform limit (after looping a few times)
if self.countonly and self.limit and self.isocount >= self.limit: return
recursion
list of donors gets smaller
list(a) constructs new list, or can use copy.copy(a)
copy doesn't need to know it's a list
deep copy is for whole data structures
don't need base case, not trying to create a single value like other recursions
builds onto a list called self.isoforms
instead when no more donors, enumerate([]) returns nothing
'''



print('#####')

seq = 'CCCCCCGTCCCCAGCCCGTCCCGTCCCAGCCCCCC'

maxs = 3
minin = 3
minex = 3
flank = 3

dons, accs = isoform2.gtag_sites(seq, flank, minex)

# test this in function and see if it returns/breaks out
enumerate([])

# recursion

def a():
    return "hello " + b()

def b():
    return "my " + c()

def c():
    return "friends"

print(a())

input = 'snake'
output = 'ekans'

def revS(s):

    if s == "":
        return ""

    return revS(s[1:]) + s[0]

print(revS('snake'))

def isPal(s):

    if len(s) == 0 or len(s) == 1:
        return True

    if s[0] == s[-1]:
        return isPal(f'{s[1:-1]}')

    return False

print(isPal('kayak'))

def d2b(dec, res):

    if dec == 0: 
        return res

    res = f'{dec%2}' + res
    return d2b(int(dec/2), res)


print(d2b(233, ''))

#30 minutes


dons = [6, 18, 24]
accs = [12, 30]

isos = []

# ians' code
def buildIsos(dons, accs, introns):

    don = dons[0]
    for aix, acc in enumerate(accs):
        if acc - don + 1 < 3: continue
        intron = (don, acc)
        iso = copy.copy(introns)
        iso.append(intron)
        isos.append(iso)

        for dix, ndon in enumerate(dons):
            elen = ndon - acc - 1
            if elen >= 3:
                ext = copy.copy(iso)
                buildIsos(dons[dix:], accs[aix:], ext)

buildIsos(dons, accs, [])
print(isos)

# sum of natural numbers
# 10: 1 + 2 + 3...+ 9 + 10 == 55

def recSum(n):

    if n <= 1:
        return n
    # shrink down the problem space
    return n + recSum(n - 1)

print(recSum(10))

# find 10 in sorted list
A = [-1, 0, 1, 2, 3, 4, 7, 9, 10, 20]

def binarySearch(A, left, right, x):

    if left > right:
        return -1
    
    mid = int((left + right) / 2)

    if x == A[mid]:
        return mid
    
    if x < A[mid]:
        return binarySearch(A, left, mid - 1, x)

    return binarySearch(A, mid + 1, right, x)

print(binarySearch(A, 0, len(A) - 1, 10))

# fibonacci 
# 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144...

def fib(n):

    if n == 0 or n == 1:
        return n
    
    else:
        return fib(n - 1) + fib(n - 2)
    
print(fib(10))

# trees 1:29:00

# binary search tree 
# print all leaf nodes
# depth first search
# optimize with memoizing and caching
# tail-call recursion