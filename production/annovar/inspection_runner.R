# -------------------------------
# THIS MAIN
# -------------------------------

offBaseThreshold <- 1

# -------------------------------
# EXPOSED API
# -------------------------------

annotatePossibleDuplicates <- function(variantsTibble) {
  tempIds <- seq(1:nrow(variantsTibble))
  variantsTibbleClone <- variantsTibble %>% useCoordinateColumns()
  variantsTibbleClone$tempId <- tempIds
  possibleDuplicatesString <- tempIds %>% lapply(function(currId) {
    currVariantTibble <- variantsTibble[currId,]
    currCoordinate <- currVariantTibble$coordinates
    currStart <- extractCoordinateStart(currCoordinate) %>% as.numeric()
    currEnd <- extractCoordinateEnd(currCoordinate) %>% as.numeric()
    possibleDuplicates <- variantsTibbleClone %>% 
      filter((abs(currStart - as.numeric(Start)) < offBaseThreshold) | (abs(currEnd - as.numeric(End)) < offBaseThreshold)) %>% 
      filter(tempId != currId)
    possibleDuplicates$tempId %>% as.character() %>% paste0(sep = ",", collapse = "")
  })
  # variantsTibbleClone <- variantsTibbleClone %>% useCoordinateStrings()
  variantsTibble$possibleDuplicates <- possibleDuplicatesString %>% unlist()
  # variantsTibbleClone %>% dplyr::select(-tempId)
  variantsTibble
}

