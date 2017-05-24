# -------------------------------
# HELPER FUNCTIONS
# -------------------------------

prepareAnnovarInputs <- function(variantsInput) {
  variantsInput %>% 
    dplyr::select(id, coordinates) %>% 
    filter(nchar(coordinates) > 0) %>%
    group_by(coordinates) %>% 
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    useCoordinateColumns() %>% 
    dplyr::select(Chr, Start, End, Ref, Alt) 
}

extractCaddPhred <- function(caddAndFred) {
  gsub("^.*,", "", caddAndFred)
}

callAnnovar <- function(tibble) {
  annovarInputFilePath <- paste0(resultsDir, geneName, "_query_annovar_input")
  write.table(tibble, file = annovarInputFilePath, quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
  annovarCommand <- paste0("perl /space/bin/annovar/table_annovar.pl ", 
                           paste0(homeDir, annovarInputFilePath), 
                           " /space/bin/annovar/humandb/ -buildver hg19 -out ", 
                           paste0(homeDir, resultsDir), geneName, 
                           "_query_annovar_output -remove -protocol refGene,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp135,ljb_all,exac03,cadd -otherinfo -operation g,r,f,f,f,f,f,f -nastring . -csvout")
  run.remote(cmd = annovarCommand, 
             remote = "apu")
  # Run annovar a second time for Phred-scaled CADD scores
  annovarCommandCaddPhred <- paste0("perl /space/bin/annovar/annotate_variation.pl ", 
                                    paste0(homeDir, annovarInputFilePath),
                                    " /space/bin/annovar/humandb -filter -dbtype cadd -buildver hg19 -out ",
                                    paste0(homeDir, resultsDir), geneName, "_query_annovar_output_cadd -otherinfo")
  run.remote(cmd = annovarCommandCaddPhred, 
             remote = "-q apu", 
             stderr.redirect = F)
  annovarResult <- read.csv(paste0(resultsDir, geneName, "_query_annovar_output.hg19_multianno.csv"))
  annovarCaddResult <- read.table(paste0(resultsDir, geneName, "_query_annovar_output_cadd.hg19_cadd_dropped"), sep="\t", stringsAsFactors = FALSE) %>% 
    dplyr::select(Chr = V3, Start = V4, End = V5, Ref = V6, Alt = V7, CADDandPhred = V2) %>% 
    mutate(CADD.phred = extractCaddPhred(CADDandPhred)) %>% 
    dplyr::select(-CADDandPhred)
  annovarResultJoined <- annovarResult %>% left_join(annovarCaddResult, by = c("Chr", "Start", "End", "Ref", "Alt"))
  annovarResultJoined
}

# -------------------------------
# EXPOSED API
# -------------------------------

annotateAnnovar <- function(variantsTibble) {
  annovarResults <- variantsTibble %>% 
    prepareAnnovarInputs() %>% 
    callAnnovar() %>%
    as_tibble()
  annovarResults$cdna <- annovarResults$AAChange.refGene %>% as.character() %>% lapply(extractCdna) %>% unlist()
  annovarResults$aachange <- annovarResults$AAChange.refGene %>% as.character() %>% lapply(extractAAChange) %>% unlist()
  annovarResults <- annovarResults %>% 
    useCoordinateStrings() %>% 
    dplyr::select(coordinates, aachange, cdna, Func.refGene, ExonicFunc.refGene, CADD.phred)
  variantsTibble %>% 
    left_join(annovarResults, by = c("coordinates"))
}



