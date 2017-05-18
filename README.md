# SFARI Variant Prioritisation Pipeline

## Introduction
Computational pipeline to identify variants of interest for functional prioritisation in SFARI collaboration's model organism systems.

## Install 

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


## Usage - setup

### SSH keys

Are you set up for passwordless ssh with keys? If no:
http://www.chibi.ubc.ca/faculty/pavlidis/wiki/display/PavLab/Connecting+to+the+Servers

Also for fixing the hostname verification failed problem that seems to crop up with Chalmers:
https://askubuntu.com/questions/45679/ssh-connection-problem-with-host-key-verification-failed-error

### Setup Config File

Edit config.R with HUGO gene symbol, specify transcript, chromosome number

### Setup Query Variants File

1. Copy template query variants file (inputs/DYRK1A_query_variants.csv)

2. Edit query variants file 

Variants file should be a .csv with nine columns:
querytype column specifies the type of query the script makes:
if "cdna":
report variant in cdna column as form 
| chr | pos | ref | alt | cdna | aachange | Source | X | querytype | denovo |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|21|NA|A860T| NA| NA| NA| NA| NA| NA|NA |


3. Save file with name GENENAME_query_variants.csv

## Usage - Variant Prioritisation

Run R scripts
```
R main.R
```

## Project Details

Main scripts:

* main.R - wrapper script - runs other scripts and sources dependencies
* config.R - loads options - edit this file for each run
* func.R - loads necessary functions
* load.R - loads necessary files - annovar, snap2, exac, etc 
* clean.R - data cleaning, harmonizing steps
* do.R - does the variant prioritisation

main project directories:

* ./ - home for the scripts and subdirectories
* ./inputs/  - general input files and query files go here
* ./outputs/  - plots + tables from variant prioritisation output here
* ./temp/ - store intermediate files - should be able to completely delete contents at any time 
