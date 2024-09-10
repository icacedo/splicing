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