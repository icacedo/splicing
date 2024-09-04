import math

# 600 at 100 simulations 0.778258 seconds
# 600 at 50 simulations 0.334279 seconds
# 600 at 200 simulations 1.321207 seconds
# isorandom gives the total number of seconds, not the average for each

# a and b values only work for dataset with 1000 tested sequences
def get_y(x, a, b):

    y = a * math.exp(b * x)

    return y/1000

a = 6.176734
b = 0.007036078

# change these variables
# l for max length to test
# n for number of simulations
l = 10000
n = 500

for i in range(500, l+50, 50):
    total_time = 0
    for j in range(300, i+10, 50):
        time = get_y(j, a, b) * n
        total_time += time
        #print(x, time)
    # i max length tested
    print(i, "\t", ((total_time/60)/60)/24)

# each time you increase the length by 100
# time to finish doubles
# hence 2 ^ n
