import subprocess

while True:
    subprocess.run('./optiso outfigs/17917.config.json --cpu 15 --program ./isoformer', shell=True)
    #sys.exit()
