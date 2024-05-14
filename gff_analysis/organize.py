import argparse
import os
import json

parser = argparse.ArgumentParser()
parser.add_argument('short_gffs', type=str, metavar='<directory>',
    help = 'directory with short abcgen gff files')
parser.add_argument('sorted_isos', type=str, metavar='<directory>',
    help = 'directory with sorted gene json files')
parser.add_argument('weights', type=str, metavar='<file>',
    help = 'file with individual gene weights and fitness')

args = parser.parse_args()

wbgs = {}
for file in os.listdir(args.short_gffs):
    gid = 'ch.'+file.split('.')[1]
    with open(args.short_gffs+file, 'r') as fp:
        for line in fp.readlines():
            line = line.rstrip()
            if 'wb id:' in line:
                wbg = line.split(' ')[3]
                wbgs[gid] = wbg


wb_isos = {}
for file in os.listdir(args.sorted_isos):
    gid = 'ch.' + file.split('.')[0]
    wb_isos[gid] = None
    with open(args.sorted_isos+file, 'r') as j:
        info = json.load(j)
        for iso in info:
            gid = iso.split('-')[0]
            if 'wb' in iso: continue
            if info[iso]['wb_frame'] == True:
                wbi = iso.split('-')[1]
                wb_isos[gid] = wbi

print(f'gid,wdpwm,wapwm,wemm,wimm,welen,wilen,icost,fitness,wbgene,wbiso')

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
            f'{gid},{wdpwm},{wapwm},{wemm},{wimm},{welen},'
            f'{wilen},{icost},{fit},{wbgs[gid]},{wb_isos[gid]}'
        )

