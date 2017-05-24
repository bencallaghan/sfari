# -------------------------------
# HELPER FUNCTIONS
# -------------------------------

loadSnap2Data <- function() {
  snap2map <- read.table(paste0(inputsDir, "snap2_fastamap"), stringsAsFactors = FALSE)
  snapfile <- snap2map[which(snap2map$V2 == as.character(geneName)), 1]
  snap2path <- paste0('/misc/pipeline42/ppdatabases/snap2results/UP000005640_9606/', snapfile, '.snap2.parsed')
  read.table(snap2path, stringsAsFactors = FALSE) %>% 
    as_tibble() %>% 
    dplyr::select(aachange = V1, snap2 = V2)
}

# -------------------------------
# THIS MAIN
# -------------------------------

snap2DataTable <- loadSnap2Data()

# -------------------------------
# EXPOSED API
# -------------------------------

annotateSnap2 <- function(variantsTibble) {
  variantsTibble %>% left_join(snap2DataTable, by = c("aachange"))
}
  
  
