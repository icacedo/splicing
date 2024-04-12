import argparse
import mdist_lib as mdl

parser = argparse.ArgumentParser()
parser.add_argument('gff1', type=str, metavar='<file>',
    help='apc or bli gff file')
parser.add_argument('gff2', type=str, metavar='<file>',
    help='wb gff file')

args = parser.parse_args()

introns1 = mdl.get_gff_intron_probs(args.gff1)
introns2 = mdl.get_gff_intron_probs(args.gff2)

print(introns1)
print(introns2)

mdist = mdl.get_mdist(introns1, introns2)

print(mdist)


