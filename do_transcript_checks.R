library(dplyr)
library(tidyr)
library(reshape2)
annotated.variants <- read.csv("temp/DYRK1A.avout.hg19_multianno.csv", header = TRUE)
annotated.variants %>% mutate(genome.change = paste0(Chr,":", Start,"-",  End,":", Ref,"/", Alt)) -> annotated.variants

# %>% filter(!duplicated(genome.change)) -> annotated.variants

annotated.variants$row <- 1:nrow(annotated.variants)

annotated.variants %>% 
  mutate(AAChange.refGene = strsplit(as.character(AAChange.refGene), ",")) %>% 
  unnest(AAChange.refGene) %>%
  mutate(AAChange.refGene = strsplit(as.character(AAChange.refGene), ";")) %>%
  unnest(AAChange.refGene) %>% 
  mutate(transcript = gsub(":.+", "", gsub("DYRK1A:","", AAChange.refGene))) %>%
  mutate(cdna = gsub(":.+", "", gsub("DYRK1A:NM_[0-9]+:.+:c.","", AAChange.refGene))) %>% 
  mutate(protein = gsub(":.+", "", gsub(".+p.","", AAChange.refGene))) %>% 
  # mutate(genome.change = paste0(Chr,":", Start,"-",  End,":", Ref,"/", Alt)) %>%
  select(c(genome.change, transcript, cdna)) -> annotated.variants.m
  

annotated.variants.m$cdna <- as.character(annotated.variants.m$cdna)
annotated.variants.m$transcript <- as.character(annotated.variants.m$transcript)

# dcast(data = annotated.variants.m, formula = annotated.variants.m$genome.change ~ annotated.variants.m$transcript, value.var = 'cdna')


k = annotated.variants.m
k.col = unique(k$transcript)
k.row = unique(k$genome.change)
k.mat = matrix('', nrow = length(k.row), ncol = length(k.col))
for (i in 1:length(k.row)) {
  subset.logical = (k.row[i] == k$genome.change)
  k.mat[i, ] = k$cdna[subset.logical][match(k.col, k$transcript[subset.logical])]
}
dimnames(k.mat) = list(k.row, k.col)

k.mat

write.csv(k.mat, "outputs/DYRK1Atranscripts_cdna", quote = FALSE)



# 
# 
# 
# annotated.variants.m %>% filter(!transcript == ".") -> annotated.variants.m
# spread(annotated.variants.m, key = transcript,value =  cdna)
# 
# # spread(annotated.variants.m,transcript, cdna)
# # From the source:
# # "subject" and "sex" are columns we want to keep the same
# # "condition" is the column that contains the names of the new column to put things in
# # "measurement" holds the measurements
# library(reshape2)
# 
# 
# data_with_index <- ddply(annotated.variants.m, .(genome.change), mutate, 
#                          index = paste0('medication', 1:length(Name)))    
# dcast(data_with_index, Name~ index, value.var = 'MedName')
# 
# 
# dcast(annotated.variants.m, row ~ transcript, value.var="cdna")
# data_wide
# 
# 
# 
# 
# iris.df = as.data.frame(iris)
# iris.df$row <- 1:nrow(iris.df)
# 
# long.iris.df = iris.df %>% gather(key = feature.measure, value = size, -Species)
# w.iris.df = long.iris.df %>% spread(key = feature.measure, value = size, -Species)
# 
# 
# 
# iris.df = as.data.frame(iris)
# iris.df$row <- 1:nrow(iris.df) ## added explicit row numbers
# long <- gather(iris.df, vars, val, -Species, -row) ## and made sure to keep em
# wide <- spread(long, vars, val)
# 



