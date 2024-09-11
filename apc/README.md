# Running the APC algorithm on our curated dataset
Download the APC dataset from github, into the appropriate location
```
├── Code/
│   ├── splicing/
│   ├── isoforms/
│   │   ├── apc/
|   |   |   ├── fastas & gffs

git clone https://github.com/KorfLab/isoforms.git
```
Create feature models from APC gene set  
*small differences between isoformer/ and splicing/ length models
Make sure you run this in a conda environment that has openturns installed
```
python3 make_models.py ../../isoforms/apc/ --outdir models/
```
Just use isoforms/models/ for now...  
Write file of commands to use for APC algorithm
```
python3 write_apc_cmds.py ../../isoforms/apc/ --weights ../../../data/1045weights.txt --outfile apc_cmds.txt --gff_out ../../../data/apc_gffs/ --gff_name 1045.apc
```
Run APC on multiple threads
```
python3 multi_apc.py apc_cmds.txt --cpus 15
```

This is how long it took spifire to run the APC algorithm on the 1045 dataset using 32 threads:  
time: 494.5164587497711  
real    8m14.742s  
user    198m58.396s  
sys     3m32.204s

This is how long it took lightning to run the APC algorithm on the 1045 dataset using 15 threads:  
time: 1047.4677140712738  
real    17m27.511s  
user    225m42.328s  
sys     1m31.272s