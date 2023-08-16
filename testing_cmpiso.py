import isoform_fixed as isof
import sys

gff = sys.argv[1]

introns = isof.get_introns(gff)
print(introns)







