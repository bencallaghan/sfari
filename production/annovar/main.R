library(tidyverse)
library(dplyr)

# -------------------------------
# CONFIGS
# -------------------------------

homeDir <- "~/ws/sfari/production/annovar/"
inputsDir <- "inputs/"
# IMPORTANT! --- change with query ---
queryDir <- paste0(inputsDir, "dyrk1a/")
# IMPORTANT! --- change with query ---
resultsDir <- paste0("results/dyrk1a_", format(Sys.time(), "%b_%d_%H:%M:%S") %>% as.character(), "/")

# -------------------------------
# LOAD BASIC PARAMETERS
# -------------------------------

source("coordinates_utils.R")
setwd(homeDir)
source(paste0(queryDir, "query_params.R"))
dir.create(resultsDir)

# -------------------------------
# LOAD RUNNERS
# -------------------------------

source("annovar_runner.R")
source("snap2_runner.R")
source("inspection_runner.R")

# -------------------------------
# MAIN
# -------------------------------

queryVariantsInput <- read_csv(paste0(queryDir, "query_variants.csv"))

annotatedVariants <- queryVariantsInput %>% annotateAnnovar() 
annotatedVariants <- annotatedVariants %>% annotateSnap2()

annotatedVariants <- annotatedVariants %>%
  dplyr::select(coordinates, aachange, cdna, 
                mutation_type, denovo, Func.refGene, ExonicFunc.refGene,
                CADD.phred, snap2, everything()) %>% 
  arrange(source_type, Func.refGene, ExonicFunc.refGene, coordinates, denovo, CADD.phred, snap2)

annotatedVariants <- annotatedVariants %>% annotatePossibleDuplicates()

# output it in csv!
annotatedVariants %>% write_csv(paste0(resultsDir, geneName, "_ANNOTATED_VARIANTS.csv"))
