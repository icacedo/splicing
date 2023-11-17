import sys

mdists = sys.argv[1]

gIDs = []
mds = []
with open(mdists, 'r') as fp:
	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('G'):
			gID = line.split(' ')[2]
			gIDs.append(gID)
		if line.startswith('m'):
			mdist = line.split('=')[1]
			mds.append(mdist)

mIDs = {}
for g, m in zip(gIDs, mds):
	mIDs[g] = m

#for i in mIDs:
#	print(mIDs[i], i)

minmd = min(mIDs.values())
print(minmd)

for m in mIDs:
	if mIDs[m] == minmd:
		print(m, mIDs[m]) 



