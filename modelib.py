import sys
import gzip

# unit testing in python
# https://www.dataquest.io/blog/unit-tests-python/

# for testing
fp = sys.argv[1]

# definitions of rectangular and triangular smoothing found here:
# https://terpconnect.umd.edu/~toh/spectrum/Smoothing.html

def read_intfile(fp):

	if fp.endswith(".gz"):
		with gzip.open(fp, 'r') as intfile:            
			for line in fp.readlines():
				line = line.rstrip()
				if isinstance(line, bytes):
					line = line.decode()
				yield line

	else:
		with open(fp, 'r') as intfile:
			for line in intfile.readlines():
				line = line.rstrip()
				yield line

# returns a histogram for frequencies and counts
def get_intbins(fp, nbins=None, prec=3):
	
	intlines = []
	total_obs= 0
	for l in read_intfile(fp):
		intlines.append(l)
		total_obs += 1

	intsizes = []
	for intron in intlines:
		intsizes.append(len(intron))
	
	intfreqs = []
	total = len(intsizes)
	for count in intsizes:
		fq = count/total
		fq2 = f"{fq:.{prec}f}"
		intfreqs.append(float(fq2))
	
	intcount_bins = [0 for x in range(max(intsizes)+1)]
	for i in range(len(intsizes)):
		intcount_bins[intsizes[i]] += 1

	intfreq_bins = []
	for i in range(len(intcount_bins)):
		freq = intcount_bins[i]/total_obs
		f = f"{freq:.{prec}f}"
		intfreq_bins.append(float(f))
	
	# only append up to nbins, for testing
	return intcount_bins[:nbins], intfreq_bins[:nbins], intsizes, intfreqs

#fbins = get_intbins(fp)[1]
#print(fbins)

# need to finish putting the smoothing code in here
# forgot i needed only a list of intron lengths for curve fitting?







