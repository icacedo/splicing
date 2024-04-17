import sys
import mdist_lib as mdl 
import isoform as iso

wb_gff = sys.argv[1]
apc_gff = sys.argv[2]


def read(gff):

	with open(gff, 'r') as fp:
		for line in fp.readlines():
			line = line.rstrip()
			line = line.split('\t')
			yield line

for line in read(wb_gff):
	if line[2] == 'intron':
		print(line)

total = 0
introns = {}
for line in read(apc_gff):
	if len(line) == 9:
		if line[2] == 'intron':
			intron = (int(line[3]), int(line[4]))
			if intron not in introns:
				introns[intron] = 0
			score = float(line[5])
			introns[intron] += score
			total += score
			#print(score)
			#print(line[2:6])

for i in introns:
	introns[i] = introns[i]/total

print(total)
print(introns)

introns1 = mdl.get_gff_intron_probs(apc_gff)
introns01 = iso.get_introns(apc_gff)

print(introns1)

print(introns01)

tot = 0
for i in introns1:
	tot += introns1[i]


print(tot)


introns2 = mdl.get_gff_intron_probs(wb_gff)
print(introns2)

mdist = mdl.get_mdist(introns1, introns2)
print(mdist)
