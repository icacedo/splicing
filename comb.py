import itertools 

p = [0, 1, 2]

for i in itertools.combinations(p, 2):
    print(i)

print('#####')

p2 = [0, 1, 2, 3]

for i in itertools.combinations(p2,3):
    print(i)

print('#####')
def comb(p, r):
    print(p[0], p[1])
    print(p[0], p[2])
    print(p[1], p[2])
    
comb(p, 2)

print('#####')

print(p2[0], p2[1], p2[2])
print(p2[0], p2[1], p2[3])

pairs = []
r = 3
for i in range(len(p2)-r+1):
    #for j in range(
    print(i)

#    for j in range(1+i, len(p2)):
 #       pair = (p2[i], p2[j])
  #      pairs.append(pair)

print('#####')

pairs = []
for i in range(len(p)):
    for j in range(1+i, len(p)):
        #print(p[i], p[j])
        pairs.append([p[i], p[j]])

print(pairs)

# this works, but no r 
for p in pairs:
    for i in range(p[-1]+1, len(p2)):
       print(p, i)
    
