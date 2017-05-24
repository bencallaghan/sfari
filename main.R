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
  dir.inputs <- paste0(dir.home, "inputs/")
  dir.outputs.parent <- paste0(dir.home, "outputs/")
  dir.outputs <- paste0(paste0(dir.home, "outputs/"), gene.i$name,"_",format(Sys.time(), '%m_%d_%H.%M'),"/")
  dir.plots <- paste0(dir.home, "plots/")
  dir.cor.plots <- paste0(dir.home, "plots/corr")
  dir.temp <- paste0(dir.home, "temp/")
} else{
  #setwd("/home/bcallaghan/OttoSu/Projects/SFARI")
  dir.inputs <- "/home/echu113/ws/sfari/inputs/"
  dir.outputs.parent <- ("/home/echu113/ws/sfari/outputs/")
  dir.outputs <- paste0("/home/echu113/ws/sfari/outputs/",gene.i$name,"_",format(Sys.time(), '%m_%d_%H.%M'),"/")
  dir.plots <- "/home/echu113/ws/sfari/plots/"
  dir.cor.plots <- "/home/echu113/ws/sfari/plots/corr"
  dir.temp <- "/home/echu113/ws/sfari/temp/"
}

path.fasta <- paste0(dir.inputs,gene.i$name,".fa")
path.annovar.out <- paste0(dir.temp,gene.i$name,"_anno.hg19_multianno.csv")
path.annovar.out2 <- paste0(dir.temp,gene.i$name,"_anno2.hg19_cadd_dropped")
path.annovar.additional.out <- paste0(dir.temp,gene.i$name,"_additional_anno.hg19_multianno.csv")
path.annovar.additional.out2 <- paste0(dir.temp,gene.i$name,"_additional_anno2.hg19_cadd_dropped")
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



