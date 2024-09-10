# Use a genetic algorithm to train feature model weights for APC
Download optiso and isoformer from github  
isoformer needs to be compiled  
Crete bin/ and lib/ directories to make isoformer and optiso executable
Add bin/ to path by adding this line to your bash.rc file
```
export PATH=$PATH:$HOME/Code/bin
export PYTHONPATH=$PYTHONPATH:$HOME/Code/lib
```
Follow the github instructions to compile isoformer (use make)
```
├── Code/
│   ├── bin/
│   |   ├── isoformer
│   |   ├── optiso
│   ├── splicing/
│   ├── isoforms/
│   │   ├── optiso
|   ├── genomikon/
│   │   ├── isoformer/

git clone https://github.com/KorfLab/isoforms.git
git clone https://github.com/KorfLab/genomikon.git

cd Code/bin/
ln -s ../genomikon/isoformer/isoformer
ln -s ../isoforms/optiso
```

isoformer and optiso should be executable from anywhere now

Note: APC gene set as of May 13, 2024 in isoformer/apc/ was used to calculate weights found in /share/korflab/home/ismael/data/optiso (1045 build)

### Manifest
optcon.json: config file with two genes for testing
optper.pl: simple script that runs optiso 
run_optiso.py: script that runs optiso
write_config.py: creates configuration files

