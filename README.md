# Modeling the stochastic nature of alternative splicing
Imperfect selection of intron splice signals by the alternative splicing machinery results in a large variety of low-abundance, erroneus splice isoforms alongside the canonical isoform. This project aims to predict the stocahstic background of alternative splicing using probabilistic models. A small subset of genes from the _C. elegans_ genome (see KorfLab/datacore/project_splicing/ to build this dataset) are used as a test set. All possible combinations of splice isoforms are created and scored using the apc algorithm.
## Organization
### ```apc/```
```apc_isogen.py```
+ runs apc algorithm on a single gene
```apc_model_lib.py```
+ contains functions that are used to generate apc isoforms and probabilistic models
```make_models.py```
+ creates .tsv files for acceptor/donor pwms, exon/intron Markov models, and exon/intron length models
```multi_apc.py```
+ parallelizes apc_isogen.py to be used on every gene in the apc dataset
```write_apc_cmds.py```
+ writes a text file with commands to be run in parallel  	   
### ```icost/```
```apc_pickler.py```
+ creates .pkl files from the output of aml.apc()
+ takes one gene as input
```apc_score.py```
+ uses single a .pkl and fasta file
+ scores each isoform and writes to stout in gff format
```icost_scoring.py```
+ returns .json file with all scored isoforms and tested icost values
+ uses apc_score.py to score isoforms
```mdist_lib.py```
+ functions to be called in icost_scoring.py
```run_apc_pickler.py```
+ creates .pkl files for all genes
+ uses apc_pickler.py
```avg_mdist.py```
+ takes .json file from icost_scoring.py
+ returns average Manhattan distance for each tested icost across all genes/isoforms
+ view output to pick best icost
### ```gff_analysis/```
code that parses through apc results and organizes isoforms
+ work in progress
### ```mkmdls_out/```
output of make_models.py in apc/, includes .tsv files for position weight matrices, Markov models and length models
### ```results/```
code that generates svg/html files to view genes and isoforms, written by Ian
### ```data/```
contains curated dataset for apc analysis, apc generated gff files and training sequences for apc models
### ```other/```
stuff im trying to learn
### ```arch/```
+ old apc code
### ```arch2/```
+ testing and devolpment scripts for rewriting apc code
## how to use code
### unzipping files to be used
```
cd data/
mkdir build/
tar -xvf file.tar.gz -C build/
gunzip file.txt.gz -c > build/file.txt
```
### train probabilistic models
```
mkdir mkmdls_out/
cd apc/
python3 make_models.py --extxt ../data/build/exon.txt --intxt ../data/build/intron.txt --len_limit 500 --dntxt ../data/build/donor.txt --actxt ../data/build/acceptor.txt --outdir ../mkmdls_out/
```
### calculate icost
```
cd icost/
ln -s ../apc/apc_model_lib.py
python3 run_apc_pickler.py ../data/build/apc/ --outdir /home/ismael/Data/
python3 icost_scoring.py /home/ismael/Data/apc_pickles/ ../data/build/apc/ --exon_len ../mkmdls_out/exon_len.tsv --intron_len ../mkmdls_out/intron_len.tsv --intron_mm ../mkmdls_out/intron_mm.tsv --exon_mm ../mkmdls_out/exon_mm.tsv --intron_mm ../mkmdls_out/intron_mm.tsv --donor_pwm ../mkmdls_out/donor_pwm.tsv --acceptor_pwm ../mkmdls_out/acceptor_pwm.tsv --icost_range_up 50 --icost_step 1
```
### generate apc isoforms
```
cd apc/
python3 write_apc_cmds.py ../data/build/apc/ --outfile apc_cmds.txt --gff_out /home/ismael/Data/apcgen_gffs/ --gff_name apcgen.ws290 --exon_len ../mkmdls_out/exon_len.tsv --intron_len ../mkmdls_out/intron_len.tsv --exon_mm ../mkmdls_out/exon_mm.tsv --intron_mm ../mkmdls_out/intron_mm.tsv --donor_pwm ../mkmdls_out/donor_pwm.tsv --acceptor_pwm ../mkmdls_out/acceptor_pwm.tsv
python3 multi_apc.py apc_cmds.txt --cpus 14
```
write_apc_cmps.py will output a text file with commands to run in multi_apc.py
make sure openturns is importable for apc_model_lib.py so apc runs correctly






