### Introduction
This software package contains the scripts that were used to analyze resistance data from different sources and build a summary table
#### The raw data (resistance, critical concentrations and strain identification)
The raw data are available on the directory ./sources
### How to generate the summary table of the resistance data
You should run the scripts in each directory
```
# Coll et al.
sbatch -p short --mem=10G --job-name="coll" -n 1 -o out.txt -t 0-12:00 ./Coll-et-al_2018_Nature-Genetics_s41588-017-0029-0/slurm_script.sh
# Internal isolates (Pools, CETR, TDR)
./internal/internal.py 
# Patric
./patric/patric.py
# Curated from the literature
sbatch -p short --mem=10G --job-name="curated" -n 1 -o out2.txt -t 0-12:00 ./curated_phenotypes/slurm_script.sh
# ReseqTB
./reseqtb/reseqtb.py
```
The scripts will generate the .tsv (intermediate table) and .res files (resistance data).





