import os

lcounts = {}
for file in os.listdir('../data/build/apc282/'):
    if file.endswith('.fa'):
        linec = 0
        name = file.split('.')[1]
        with open(f'../data/build/apc282/{file}', 'r') as fp:
            for fp in fp.readlines():
                linec += 1
        lcounts[name] = linec


for i in lcounts:
    print(i, lcounts[i])

print(min(lcounts.values()))

shorts = [key for key in lcounts if lcounts[key] == 6]

print(shorts)