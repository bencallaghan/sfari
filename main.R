# SFARI gene and variant prioritisation
# main.R is setup and source calls for all associated scripts
# func.R and do.R have more detail for individual functions/tasks

cat("###Starting Variant Prioritisation###\n") # Progress message

# Setup: Dependencies -----------------------------------------------------

# source("http://bioconductor.org/biocLite.R")
# library(BiocInstaller)
# biocLite("S4Vectors")
# biocLite("Biostrings")
# biocLite('biomaRt')

# install.packages('stringr')
# install.packages('seqinr')
# install.packages('ssh.utils')
# install.packages('biomaRt')
# install.packages('dplyr')
# install.packages('xtable')

library('Biostrings')
library(stringr)
library(seqinr)
library(ssh.utils)
library('biomaRt')
library(dplyr)
library(xtable)
library(gridExtra)
library(ggplot2)
library(knitr)
library(GGally)


# Setup: Options ----------------------------------------------------------

# Gene is either a single string or a vector of strings (a gene list)
# gene.list <- "PTEN" 

# Setup: Variant Prioritisation run specifications

source("config.R")


# Setup: Environment ------------------------------------------------------
# Should setwd() in project directory
# Need to change directories to match your own directory structure
if(opt.session.local == FALSE){
  #setwd("/home/bcallaghan/Projects/SFARI/")
  dir.home <- "/home/bcallaghan/Projects/SFARI/"
  dir.inputs <- "/home/bcallaghan/Projects/SFARI/inputs/"
  dir.outputs.parent <- ("/home/bcallaghan/Projects/SFARI/outputs/")
  dir.outputs <- paste0("/home/bcallaghan/Projects/SFARI/outputs/",gene.i$name,"_",format(Sys.time(), '%m_%d_%H.%M'),"/")
  dir.plots <- "/home/bcallaghan/Projects/SFARI/plots/"
  dir.cor.plots <- "/home/bcallaghan/Projects/SFARI/plots/corr"
  dir.temp <- "/home/bcallaghan/Projects/SFARI/temp/"
} else{
  #setwd("/home/bcallaghan/OttoSu/Projects/SFARI")
  dir.home <- "/home/bcallaghan/OttoSu/Projects/SFARI/"
  dir.inputs <- "/home/bcallaghan/OttoSu/Projects/SFARI/inputs/"
  dir.outputs.parent <- ("/home/bcallaghan/OttoSu/Projects/SFARI/outputs/")
  dir.outputs <- paste0("/home/bcallaghan/Projects/SFARI/outputs/",gene.i$name,"_",format(Sys.time(), '%m_%d_%H.%M'),"/")
  dir.plots <- "/home/bcallaghan/OttoSu/Projects/SFARI/plots/"
  dir.cor.plots <- "/home/bcallaghan/OttoSu/Projects/SFARI/plots/corr"
  dir.temp <- "/home/bcallaghan/OttoSu/Projects/SFARI/temp/"
}

path.fasta <- paste0(dir.inputs,gene.i$name,".fa")
path.annovar.out <- paste0(dir.temp,gene.i$name,"_anno.hg19_multianno.csv")
path.annovar.out2 <- paste0(dir.temp,gene.i$name,"_anno2.hg19_cadd_dropped")
path.pp.out <- paste0(dir.inputs, gene.i$name, "PP.tsv")
path.marv <- paste0(dir.inputs,"MARV_ASD_muts_hg19_multianno.csv")
path.anno.in <- paste0(dir.temp,as.character(gene.i$name),"_anno_in")
path.gene.metrics <- paste0(dir.outputs.parent,"gene_level_info")
path.gene.lit.variants <- paste0(dir.inputs, gene.i$name, ".lit.variants")
path.query.variants <- paste0(dir.inputs, gene.i$name, "_query_variants.csv")
path.log <- paste0(dir.outputs,"log.txt")

dir.create(dir.outputs) # create a new directory (timestamped) for every run
log_con <- file(path.log,open="a")

# Setup: Tests ------------------------------------------------------------


# Source Calls ------------------------------------------------------------

source("func.R")
source("load.R")
source("clean.R")
source("do.R")



