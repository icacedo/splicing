import sys
import os

# directory with abc generated gff files
gdir = sys.argv[1]

if not os.path.exists('build/20isos/'): os.makedirs('build/20isos/')

for file in os.listdir(gdir):
    gid = file.split('.')[1]
    new = f'build/20isos/ch.{gid}.abcgen.282.20.gff'
    with open(gdir+file, 'r') as fp:
        if os.path.isfile(new): os.remove(new)
        count = 0 
        for line in fp.readlines():
            with open(new, 'a') as gff:
                line = line.rstrip()
                gff.write(f'{line}\n')
                if line == '': count += 1
                if count == 21: break
    
    


