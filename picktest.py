import pickle
import sys

name = sys.argv[1]

with open(name, 'rb') as pick:
	isos = pickle.load(pick)

count = 0
for i in isos:
	print(i)
	count += 1

print(count)
