# Code for modeling the stochastic background of RNA splicing...

## Questions for later ##
+ exon EVD scores exons <20 bp as a non-zero probability
+ there are no exons <20 bp in exon.txt
+ g parameter is negative, so x=0 (exon len) will be greater
+ see memoize_fdist in modelib
+ probability/score should be 0/-100?

## MANIFEST ##
### splicing/
+ modelib.py - apc function library
+ make_models.py - generates probabilistic models
+ pwm_scoring.py - scores donor and acceptor sites with pwm
+ isotyper.py - categorizes isoforms
### splicing/data/
+ copied from arch/
+ apc dataset, exon/intron sequences, donor/acceptor sequences
### splicing/mkmdls_out/
+ output of make_models.py
+ .tsv files for each model
### splicing/results/
+ output from Ian's code
+ svg files
+ kinda forget how this works
### splicing/arch/
+ isoform.py - library for most functions
+ geniso - generates isoforms
+ cmpiso - compares isoforms
+ optiso - optimizes parmeters for geniso
+ randomly.py - examines splicing in random sequences
+ data - directory training data and scripts
	+ extract_sequences.py - uses Lyman2020 favorites
	+ build_models.py - creates pwms, etc
 ### splicing/arch2/
 + devlopment scripts for apc rewrite



## To Do ##
- Compare Ian's apc code to rewrite
- Compare RNA-seq data with simulated data
- Compare RNA-seq data with gff-generated isoforms (biologically likely)
- Compare simulated data with gff-generated isoforms (biologically likely)
- Describe differences between APC 668 and WormBase, see (ask permission to view):
	https://docs.google.com/spreadsheets/d/1tzfDo145Nt0MtEyEqy__WnaV-PtKA79LN6mrsWuZbt0/edit#gid=0
