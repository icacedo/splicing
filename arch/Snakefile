# use --delete-all-output to remove all files in pipeline
# used this command to get list of ids:
# find data/apc/ -name "*.fa" | cut -d . -f 2 > id_list.txt
'''
fp = open('id_list.txt')
IDs = []
for id in fp.readlines():
	IDs.append(id.rstrip())
'''

configfile: 'config.yml'
IDs = config['IDs']

min_intron = config['min_intron']
min_exon = config['min_exon']
max_splice = config['max_splice']
flank = config['flank']

dpwm = config['dpwm_path']
apwm = config['apwm_path']
emm = config['emm_path']
imm = config['imm_path']
elen = config['elen_path']
ilen = config['ilen_path']

limit = config['iso_limit']

rule all:
	input:
		expand('cmpiso_out/ch.{id}.manhattan', id=IDs)
	
rule run_geniso:
	input:
		f1=('test_data/ch.{sample}.fa'),
		f2=('test_data/ch.{sample}.gff3')
	output:
		'geniso_out/ch.{sample}.apc.iso.gff'
	shell:
		'./geniso {input.f1} --min_intron {min_intron} --min_exon {min_exon} '
		'--max_splice {max_splice} --flank {flank} --dpwm {dpwm} ' 
		'--apwm {apwm} --emm {emm} --imm {imm} --elen {elen} --ilen {ilen} '
		'--introns {input.f2} > {output}'
		
rule run_cmpiso:
	input:
		'geniso_out/ch.{sample}.apc.iso.gff',
		'test_data/ch.{sample}.gff3'
	output:
		'cmpiso_out/ch.{sample}.manhattan'
	shell:
		'./cmpiso {input} > {output}'






