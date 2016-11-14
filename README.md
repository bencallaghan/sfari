# sfari project R scripts


R code for the SFARI project

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
