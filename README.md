# Code for modeling the stochastic background of RNA splicing...

## MANIFEST ##

+ isoform.py - library for most functions
+ geniso - generates isoforms
+ cmpiso - compares isoforms
+ optiso - optimizes parmeters for geniso
+ randomly.py - examines splicing in random sequences
+ data - directory training data and scripts
	+ extract_sequences.py - uses Lyman2020 favorites
	+ build_models.py - creates pwms, etc

## To Do ##

- Compare RNA-seq data with simulated data
- Compare RNA-seq data with gff-generated isoforms (biologically likely)
- Compare simulated data with gff-generated isoforms (biologically likely)
- Describe differences between APC 668 and WormBase, see (ask permission to view):
	https://docs.google.com/spreadsheets/d/1tzfDo145Nt0MtEyEqy__WnaV-PtKA79LN6mrsWuZbt0/edit#gid=0
