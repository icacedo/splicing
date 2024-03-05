# modeling the stochastic nature of alternative splicing

## directory structure
need apc pickler
then icost scoring
then avg mdis


### arch/
+ old apc code
### arch2/ 
+ testing and devolpment scripts for rewriting apc code
## Directions for icost testing ##
+ apc_pickler.py
  	+ single apc gene fasta file as input
	+generates apc isoform .pkl files
  	+ places .pkl file in apc_pickles/
+ apc_score.py
  	+ scores isoforms in apc .pkl files
  	+ input single .pkl file
  	+ output single .gff file
+ icost_testing.py
  	+ input directory with apc pickle files
  	+ input directory with apc fasta files
  	+ input scoring model .tsv files
  	+ imports apc_score.py
     
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
