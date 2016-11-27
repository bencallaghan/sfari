# SFARI gene and variant prioritisation
# main.R is setup and source calls for all associated scripts
# func.R and do.R have more detail for individual functions/tasks



# Setup: Dependencies -----------------------------------------------------

source("https://bioconductor.org/biocLite.R")
# biocLite('biomaRt')

install.packages('stringr')
install.packages('seqinr')
install.packages('ssh.utils')
install.packages('biomaRt')
install.packages('dplyr')
install.packages('xtable')

# 




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


gene.list <- data.frame(name = c("PTEN","TBR1", "GRIN2B", "DYRK1A","SYNGAP1"),
                         transcript = c("NM_000314","NM_006593","NM_000834","NM_001396","NM_006772"),
                         chromosome = c(10,2,12,21,6))
i <- 4
gene.i <- gene.list[i,]

opt.annovar.cache = TRUE # If TRUE, check for cached copy of annovar in temp before rerunning 
opt.generanks.cache = TRUE # If TRUE, check for existing final gene ranking
opt.session.local = FALSE # If TRUE, change all dirs accordingly

# Setup: Environment ------------------------------------------------------

if(opt.session.local == FALSE){
  setwd("/home/bcallaghan/Projects/SFARI/")
  dir.home <- "/home/bcallaghan/Projects/SFARI/"
  dir.inputs <- "/home/bcallaghan/Projects/SFARI/inputs/"
  dir.outputs <- "/home/bcallaghan/Projects/SFARI/outputs/"
  dir.plots <- "/home/bcallaghan/Projects/SFARI/plots/"
  dir.cor.plots <- "/home/bcallaghan/Projects/SFARI/plots/corr"
  dir.temp <- "/home/bcallaghan/Projects/SFARI/temp/"
} else {
  setwd("/home/bcallaghan/OttoSu/Projects/SFARI")
  dir.home <- "/home/bcallaghan/OttoSu/Projects/SFARI/"
  dir.inputs <- "/home/bcallaghan/OttoSu/Projects/SFARI/inputs/"
  dir.outputs <- "/home/bcallaghan/OttoSu/Projects/SFARI/outputs/"
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
path.gene.metrics <- paste0(dir.outputs,"gene_level_info")
path.gene.lit.variants <- paste0(dir.inputs, gene.i$name, ".lit.variants")

# Setup: Tests ------------------------------------------------------------


# Source Calls ------------------------------------------------------------

source("func.R")
source("load.R")
source("clean.R")
source("do.R")



