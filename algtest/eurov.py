import random 
import statistics

scores = [12, 10, 8, 7, 6, 5, 4, 3, 2, 1]
scores = [12, 10, 8, 7, 0, 0, 0, 0]
countries = [x for x in range(6)]
judges = 4
trials = 2

datasets = {}
for j in range(judges):
	random.shuffle(countries)
	for k in range(len(countries)):
		if k > len(scores)-1:
			score = 0
		else:
			score = scores[k]
		if countries[k] not in datasets:
			datasets[countries[k]] = []
			datasets[countries[k]].append(score)
		else:
			datasets[countries[k]].append(score)


print(datasets)



'''
sims = {}
for i in range(trials):
	all_scores = {}
	for j in range(3):
		random.shuffle(countries)
		judged = {}
		for k in range(len(countries)):
			if k > len(scores)-1:
				score = 0
			else:
				score = scores[k]
			judged[countries[k]] = score
		all_scores[j] = judged
	sims[i] = all_scores
	
datasets = {}
for sim in sims:	
	dataset = {}
	for item in sims[sim].items():
		print(item)
		for country in item[1]:
			if country not in dataset:
				dataset[country] = []
				dataset[country].append(item[1][country])
			else:
				dataset[country].append(item[1][country])
'''
