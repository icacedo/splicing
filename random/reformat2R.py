# reformat isorandom output files for R readability

import sys
import os
import csv

ran_dir = sys.argv[1]

times = {}
data = []
for file in os.listdir(ran_dir):
    id = file.split('.')[0]
    with open(f'{ran_dir}/{file}', 'r') as fp:
        for line in fp.readlines():
            line = line.rstrip()
            if line.startswith('#'):
                times[id] = line.split(' ')[3]
                continue
            data.append(line)

with open('times.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    for id in sorted(times):
        writer.writerow([id,times[id]])

with open('randata.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    for d in data:
        d = d.split('\t') 
        writer.writerow([d[0], d[1], d[2], d[3]])