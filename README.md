### Introduction
This software package contains the scripts that were used to analyze resistance data from different sources and build a summary table

#### The raw data (raw resistance data and critical concentrations)
The raw data are available on the directory ./sources.

### How to generate the summary table of the resistance data
You should run the scripts in each directory
```
# Coll et al.
sbatch -p short --mem=10G --job-name="coll" -n 1 -o ./Coll-et-al_2018_Nature-Genetics_s41588-017-0029-0/output_slurm_script.txt -t 0-12:00 ./Coll-et-al_2018_Nature-Genetics_s41588-017-0029-0/slurm_script.sh
# Internal isolates (Pools, CETR, TDR)
./internal/internal.py
# Patric
./patric/patric.py
# Curated from the literature
sbatch -p short --mem=10G --job-name="curated" -n 1 -o ./curated_phenotypes/out_slurm_script.txt -t 0-12:00 ./curated_phenotypes/slurm_script.sh
# ReseqTB
./reseqtb/reseqtb.py
```
The scripts will generate the .tsv (intermediate table -- for some of the datasets) and .res files (resistance data -- for all the datasets).
Now you just need to open a python shell and run:
```
from utils.generate_summary import *

generate_summary_from_res_remove_collisions(["./Coll-et-al_2018_Nature-Genetics_s41588-017-0029-0/Coll-et-al_2018_Nature-Genetics_s41588-017-0029-0.res","./internal/internal.res","curated_phenotypes/curated_phenotypes.res","./patric/patric.res","./reseqtb/reseqtb.res"],"summary_table_resistance-<version>.tsv")
```

NOTE: of course you need to increment the version number if you made some modifications to the scripts. The versioning works like that:
* X.Y: increment X if you add a new set of strains with resistance_data; increment Y if you improve the quality of the current data set of strains.
