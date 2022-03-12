rule all:
	input:
		"ch.7402.apc.iso"
		
rule run_geniso:
	input:
		"data/apc/ch.7402.fa"
	output:
		"ch.7402.apc.iso"
	shell:
		"./geniso data/apc/ch.7402.fa > ch.7402.apc.iso"
