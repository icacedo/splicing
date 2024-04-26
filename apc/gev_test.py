import apc_model_lib as aml
#import openturns as ot
import sys

# remove import openturns from acp_model_lib

seqs = sys.argv[1]

reseqs = aml.read_txt_seqs(seqs)

data = aml.get_exinbins(reseqs)[2]

test = [x for x in data if x < 500]

print(test)

#data = aml.fdist_params(reseqs, size_limit=500)

#y1, y2 = aml.memoize_fdist(data, a, b, g, size_limit)

#print(y1)
#rint(y2)
