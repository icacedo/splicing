# scoring program for intron length model
# graph visualization in curve_fitting.ipynb
# memoize: calculate pdf once, put the output in an array for lookup
# this code will go in modelib.py

import openturns as ot
import openturns.viewer as viewer
import math
import modelib
import sys

fp = sys.argv[1]

data = modelib.get_intbins(fp, 500, 5)[2]
print(data)
size_limit = 250

sample = ot.Sample([[x] for x in data if x < size_limit])

distFrechet = ot.FrechetFactory().buildAsFrechet(sample)

a = distFrechet.getAlpha()
b = distFrechet.getBeta()
g = distFrechet.getGamma()

def frechet_pdf(x, a, b, g):
	if x < g: return 0
	z = (x-g)/b
	term1 = (a/b)
	term2 = z**(-1-a)
	term3 = math.exp(-z**-a)
	return term1 * term2 * term3

x_values = []
y_values = []
y_scores = []
for i in range(size_limit):
	x_values.append(i)
	y = frechet_pdf(i, a, b, g)
	y_values.append(y)
	if y == 0: y_scores.append(-99)
	else: y_scores.append(math.log2(y))

# this is the score lookup dictionary
data = {
	'x': x_values,
	'y': y_values
}

ilen = 101
print(data['y'][ilen])







