import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('gdir', type=str, metavar='<director>',
    help = 'directory with short abcgen gff files')
parser.add_argument('weights', type=str, metavar='<file>',
    help = 'file with individual gene weights and fitness')

args = parser.parse_args()

'''
wbgs = {}
for file in os.listdir(args.gdir):
    gid = 'ch.'+file.split('.')[1]
    with open(args.gdir+file, 'r') as fp:
        for line in fp.readlines():
            line = line.rstrip()
            if 'wb id:' in line:
                wbg = line.split(' ')[3]
                wbgs[gid] = wbg

print(f'gid,wdpwm,wapwm,wemm,wimm,welen,wilen,icost,fitness,wbgene')
with open(args.weights, 'r') as fp:
    for line in fp.readlines():
        line = line.rstrip()
        line = line.split('\t')
        fit = line[0]
        wdpwm = line[1]
        wapwm = line[2]
        wemm = line[3]
        wimm = line[4]
        welen = line[5]
        wilen = line[6]
        icost = line[7]
        gid = line[8]
        print(
            f'{gid},{wdpwm},{wapwm},{wemm},{wimm},'
            f'{welen},{wilen},{icost},{fit},{wbgs[gid]}'
        )
'''

#add if there are matching wormbase isoforms to list of genes