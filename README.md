# Code for modeling the stochastic background of RNA splicing...

## MANIFEST ##

+ isoform.py - library for most functions
+ asg.py - alternative splice generator
+ randomly.py - examines splicing in random sequences
+ data - directory training data and scripts
	+ extract_sequences.py - uses Lyman2020 favorites
	+ build_models.py - creates pwms, etc

## To Do ##

Comparison program needs to be able to compare real and simulated expression
patterns. Since RNA-seq data is intron-based and not isoform-based, the
comparisons have to break apart all the transcripts and keep track of the
weights of each intron.

- Compare annotated genes with RNA-seq data
- Compare RNA-seq data with simulated data
- Compare simulated data to simulated data

Preliminary results show that some artificial isoforms rank higher than real
ones. This might be solved by using weights or costs

+ Set weight for each model
+ Set cost for splicing
- Genetic algorithm for finding optimal weights and cost
