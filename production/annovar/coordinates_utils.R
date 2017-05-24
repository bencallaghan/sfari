# input - tibble containing Chr, Start, End, Ref, Alt
useCoordinateStrings <- function(tibble) {
  tibble %>% 
    mutate(coordinates = paste0(Chr, ":", Start, "-", End, ":", Ref, "/", Alt)) %>% 
    dplyr::select(-Chr, -Start, -End, -Ref, -Alt) %>%
    dplyr::select(coordinates, everything())
}

# input - tibble containing coordinates in the form of strings
# outputs opposite of convertToCoordinateStrings
useCoordinateColumns <- function(tibble) {
  tibble %>% 
    mutate(Chr = extractCoordinateChrom(coordinates),
           Start = extractCoordinateStart(coordinates),
           End = extractCoordinateEnd(coordinates),
           Ref = extractCoordianteRef(coordinates),
           Alt = extractCoordianteAlt(coordinates)) %>% 
    dplyr::select(-coordinates) %>% 
    dplyr::select(Chr, Start, End, Ref, Alt, everything())
}

extractCoordinateChrom <- function(coordinateString) {
  gsub(":.*:.*$", "", coordinateString)
}

extractCoordinateStart <- function(coordinateString) {
  stringClone <- gsub("^[0-9]*:", "", coordinateString)
  stringClone <- gsub("-.*$", "", stringClone)
  stringClone
}

extractCoordinateEnd <- function(coordinateString) {
  stringClone <- gsub("^[0-9]*:[0-9]*-", "", coordinateString)
  stringClone <- gsub(":.*$", "", stringClone)
  stringClone
}

extractCoordianteRef <- function(coordinateString) {
  stringClone <- gsub("^.*:", "", coordinateString)
  stringClone <- gsub("/.*$", "", stringClone)
  stringClone
}

extractCoordianteAlt <- function(coordinateString) {
  stringClone <- gsub("^.*:", "", coordinateString)
  stringClone <- gsub("^.*/", "", stringClone)
  stringClone
}

extractCdna <- function(variantStrings) {
  transcriptRegex <- paste0("^", geneName, ":", transcript, ":.*")
  variantString <- grep(transcriptRegex, strsplit(variantStrings %>% as.character(), ",") %>% unlist(), value = TRUE)
  cdnaString <- gsub(paste0("^", geneName, ":", transcript, ":exon[0-9]*:c\\."), "", variantString)
  cdnaString <- gsub(paste0(":.*$"), "", cdnaString)
  if (length(cdnaString) > 0 && nchar(cdnaString[[1]]) > 0) {
    return(cdnaString)
  }
  return(".")
}

extractAAChange <- function(variantStrings) {
  transcriptRegex <- paste0("^", geneName, ":", transcript, ":.*")
  variantString <- grep(transcriptRegex, strsplit(variantStrings %>% as.character(), ",") %>% unlist(), value = TRUE)
  aaChangeString <- gsub(paste0("^", geneName, ":", transcript, ".*:p\\."), "", variantString)
  if (length(aaChangeString) > 0 && nchar(aaChangeString[[1]]) > 0) {
    return(aaChangeString)
  }
  return(".")
}
