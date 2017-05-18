# SFARI Variant Prioritisation Pipeline

## Introduction
Computational pipeline to identify variants of interest for functional prioritisation in SFARI collaboration's model organism systems.

# Install 

On pavlab servers, clone the project:

```
mkdir sfari/
cd sfari/
git clone https://github.com/bencallaghan/sfari
```


Fetch updated files from MARVdb:
```
sh dump_marvdb.sh
```


# Usage

### Config File
Edit config.R with HUGO gene symbol, specify transcript, chromosome number

### Query Variants File
Copy template query file (inputs/DYRK1A_query_variants.csv)
Edit variants file 
Save file with name GENENAME_query_variants.csv


Run R scripts
```
R main.R
```

5 main scripts:

* main.R - "wrapper" script - runs everything else, loads dependencies, can make loops for running > 1 gene
* func.R - contains all functions - run this before everything else to load necessary functions!
* load.R - loads necessary files - annovar, snap2, exac, etc etc
* clean.R - data cleaning, harmonizing steps
* do.R - does the meat of the processing work

main project directories:

* ./ - home for the scripts and subdirectories
* ./inputs/  - where to put gene-specific and general input files - load them with load.R
* ./outputs/  - where a lot of plots, tables, etc spits out. Can organise into gene - specific folders
* ./temp/ - store intermediate files - should be able to completely delete contents at any time with no ill effects (other than potentially having to regenerate whatever files were inside
