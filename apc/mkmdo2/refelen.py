import sys

fp = sys.argv[1]
'''
with open(fp, 'r') as model:
    for line in model.readlines():
        line = line.rstrip()
        print(f'{line}\t1.1')
'''

with open(fp, 'r') as model:
    for line in model.readlines():
        line = line.rstrip()
        if line.startswith('%'): print(line)
        else:
            line = line.split(' ')
            if line == ['']: continue
            print(f'{line[0]}\t{line[1]}\t1.1')