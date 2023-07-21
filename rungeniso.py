import subprocess

# to be run from splicing root

files = ['ch.13302.fa', 'ch.302.fa']
geniso = 'arch/geniso'
outdir = 'build/out'

for file in files:
	ff = f'build/apc/{file}'
	cli = f'{geniso} {ff}'
	
	# if you just want to run it and save output somewhere
	subprocess.run(f'{geniso} {ff} --min_exon 100 > {outdir}/{file}.geniso', shell=True)
	
	# if you want to parse live...
	#for line in subprocess.run(f'{geniso} {ff} --min_exon 100', shell=True, \
	#	capture_output=True).stdout.decode().split('\n'):
	#	if line.startswith('#'): continue
	#	f = line.split()
	#	if len(f) == 0: continue
	#	print(f)

