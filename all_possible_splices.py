import itertools

don = [50, 100, 150];
acc = [25, 90, 120];

for d in range(0, len(don)):
	for dsites in itertools.combinations(don, d):
		for a in range(0, len(acc)):
			for asites in itertools.combinations(acc, a):
			
				# same number of sites both sides required
				if len(dsites) != len(asites): continue
				
				# make introns
				introns = []
				for i in range(len(dsites)):
					if dsites[i] < asites[i]:
						introns.append((dsites[i], asites[i]))
				if len(introns) == 0: continue
				
				# check exons at some point
				
				# print introns
				print(introns);
				
