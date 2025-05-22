import random 
import statistics

scores = [12, 10, 8, 7, 6, 5, 4, 3, 2, 1]
countries = [x for x in range(26)]
judges = 37 * 2
trials = 1000

winners = []
agg_data = []
for t in range(trials):
	dataset = {}
	for j in range(judges):
		random.shuffle(countries)
		for k in range(len(countries)):
			if k > len(scores)-1:
				score = 0
			else:
				score = scores[k]
			if countries[k] not in dataset:
				dataset[countries[k]] = 0
				dataset[countries[k]] += score
			else:
				dataset[countries[k]] += score
	winners.append(max(dataset.values()))
	for tscore in dataset.values():
		agg_data.append(tscore)

print('### all countries ###')
print('mean: ', statistics.mean(agg_data))
print('stdev: ', statistics.stdev(agg_data))
print('### winners #########')
print('mean: ', statistics.mean(winners))
print('stdev: ', statistics.stdev(agg_data))






'''
trial_sums = {}
for t in range(trials):
	dataset = {}
	for j in range(judges):
		random.shuffle(countries)
		for k in range(len(countries)):
			if k > len(scores)-1:
				score = 0
			else:
				score = scores[k]
			if countries[k] not in dataset:
				dataset[countries[k]] = []
				dataset[countries[k]].append(score)
			else:
				dataset[countries[k]].append(score)
	trial_sums.append(dataset)
	

print(dataset)
averges = {}
for data in dataset.items():
	print(data)
	print(statistics.mean(data[1]))
	sum_scores = lambda my_list: sum(my_list)
	print(sum_scores[data[1]])
'''
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
