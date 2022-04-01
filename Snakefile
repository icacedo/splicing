# use --delete-all-output to remove all files in pipeline
# used this command to get list of ids:
# find data/apc/ -name "*.fa" | cut -d . -f 2 > id_list.txt
'''
fp = open('id_list.txt')
IDs = []
for id in fp.readlines():
	IDs.append(id.rstrip())
'''
#IDs = ['4884', '7395']
configfile: 'config.yml'
IDs = config['IDs']
print('Config is: ', config)
print(IDs)

rule all:
	input:
		expand('cmpiso_out/ch.{id}.manhattan', id=IDs)
	
rule run_geniso:
	input:
		'test_data/ch.{sample}.fa'
	output:
		'geniso_out/ch.{sample}.apc.iso.gff'
	shell:
		'./geniso {input} > {output}'

rule run_cmpiso:
	input:
		'geniso_out/ch.{sample}.apc.iso.gff',
		'test_data/ch.{sample}.gff3'
	output:
		'cmpiso_out/ch.{sample}.manhattan'
	shell:
		'./cmpiso {input} > {output}'






