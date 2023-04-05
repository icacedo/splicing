import sys
import gzip

fp = sys.argv[1]

def read_file(fp):

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

xlines = 5

intlines = []
count = 0
for l in read_file(fp):
	if count <= xlines:
		intlines.append(l)
		count += 1	

intlines = [
	'AAATGCTATA',
	'ATCTA',
	'CGACTCT',
	'ACTACGTACG',
	'ACTAGG',
	'CGCGATCGTT',
	'AGCGATTTACAGCG',
	'AGCAC',
	'CTACTG',
	'AGAG',
	'AGCT'
]

print(intlines)	

intsizes = []
for i in intlines:
	intsizes.append(len(i))
print(intsizes)

intbins = [0 for x in range(max(intsizes)+1)]
print(intbins)

for i in range(len(intsizes)):
	intbins[intsizes[i]] += 1

print(intbins)
		
print('***********')
	
# try different smoothing techniques
# rectangular and with a slope
# or parabolic (goes to 0)
# start with linear vs rectangular
# rec, parameter is the length n of the rectangle
# linear, parameter is the slope m
# can be done in frequencies or sizes
# slope can create fractions
# can create a slope that goes by a set distance
# slope m, would be dependent on distance n
# how much distance? don't want to say intron of 0 is possible
# resources:
# https://terpconnect.umd.edu/~toh/spectrum/Smoothing.html
# http://sepwww.stanford.edu/sep/prof/pvi/zp/paper_html/node6.html#SECTION00121000000000000000
# http://sepwww.stanford.edu/sep/prof/pvi/zp/paper_html/node7.html

# m must be odd, smooth width
m = 5
'''
intbins = [1,2,3,4,5,6,7,8,9]
for i in range(len(intbins)):
	m2 = int((m/2) + 0.5 - 1)
	print(intbins[i-m2:i])
	print(intbins[i])
	print(intbins[i+1:i+m2+1])
	print('***')
'''
smoodata = []
for i in range(len(intbins)):
	m2 = int((m/2) + 0.5 - 1)
	bef = intbins[i-m2:i]
	now = intbins[i]
	aft = intbins[i+1:i+m2+1]
	total = sum(bef) + now + sum(aft)
	smoopt = total/m
	smoodata.append(smoopt)
	print(bef, now, aft, total, smoopt)

print(smoodata)
print('***********')
# two different ways of rectangular smoothing?	
smoodata2 = []
for i in range(len(intbins)):
	m2 = int((m/2) + 0.5 - 1)
	bef = intbins[i-m2:i]
	now = intbins[i]
	aft = intbins[i+1:i+m2+1]
	total = 0
	for n in bef:
		total += n + now
	for n in aft:
		total += n + now
	total += now + now
	smoopt2 = total/m
	smoodata2.append(smoopt2)
	print(bef, now, aft, total, smoopt2)

print(smoodata2)


