import modelib as ml
import sys

fp = sys.argv[1]

intbins = ml.get_intbins(fp)[0]

print(ml.rec_smoo(intbins))

print(ml.tri_smoo(intbins))



