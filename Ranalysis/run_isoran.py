import subprocess
import multiprocessing
import os

threads = os.cpu_count() - 1

def worker(cmd):

    result = subprocess.run(cmd, shell=True, 
                            capture_output=True).stdout.decode()
    return result.rstrip()

out = "isoran1_2/"

os.makedirs(out, exist_ok=True)

jobs = []
for i in range(1, 2+1):
    for j in range(500, 2000+50, 50):
        cmd = f"isorandom {j} 1000 --max_splice {i} > {out}rnd.{i}.{j}.txt"
        jobs.append(cmd)

pool = multiprocessing.Pool(threads)
pool.map(worker, jobs)

'''
cmds = {}
for i in range(1, 2+1):
    for j in range(500, 2000+50, 50):
        cmd = f"isorandom {j} 1000 --max_splice {i}"
        cmds[f"rnd.{i}.{j}.txt"] = cmd

for c in cmds:
    result = subprocess.run(cmds[c], shell=True, 
                            capture_output=True).stdout.decode()
    print(result.rstrip())
    with open(c, 'w') as file:
        file.write(result)
    break
'''
