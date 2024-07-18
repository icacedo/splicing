import sys
import isomod as im
import itertools as it

fasta = sys.argv[1]
seq = im.read_fasta(fasta)
dons, accs = im.get_gtag(seq[1], 100, 25)

array = [1, 2, 3, 4, 5]

# following this tutorial:
# https://jarednielsen.com/algorithm-combinations/

for i in range(len(array)):
    head = array[i]
    tail = array[i+1:]
    for j in range(len(tail)):
        h2 = tail[j]
        t2 = tail[j+1:]
        for k in range(len(t2)):
            combo = [head, tail[j], t2[k]]
            print(combo)     

def combinations(n, k):
    combos = []
    #head = []
    #tail = []

    if k == 1:
        return n
    
    for i in range(len(n)):
        head = n[i]
        tail = n[i+1:]
        for j in range(len(tail)):
            combo = [head, tail[j]]
            combos.append(combo)
    
    return combos

    
print(combinations(array, 2))

def combi(n, k):
    combos = []

    if k == 1:
        return n
    
    for i in range(len(n)):
        head = n[i]
        tail = combi(n[i+1:], k-1)
        for j in range(len(tail)):
            combo = [head, tail[j]]
            combos.append(combo)
    
    return combos

print('###')

c = combi(array, 4)
print(c)


def combation(n, k):
    combos = []

    if (k == 1):
        return n 

    for i in range(len(n)): 
        head = n[i:i+1]

        tail = combinations(n[i+1:],k-1)

        for j in range(len(tail)):
            if (type(tail[j]) == int):
                combo = head + [tail[j]]
            else:
                combo = head + tail[j]
            combos.append(combo)
    
    return combos

print('###')

x = combation(array, 4)
print(x)

print('###')

g = []
for co in it.combinations(array, 4):
    g.append(co)

print(g)

import time
start = time.time()
for i in range(len(dons)):
    it.combinations(dons, i)
end = time.time()
print(end-start)

# very slow
'''
start = time.time()
for i in range(len(dons)):
    combi(dons, i)
end = time.time()
print(end-start)
'''