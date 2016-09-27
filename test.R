
# non ranscripts stuff
BED
vars.filtered


test.non.transcript <- function(BED, vars.filtered){
  print(paste("Before exon 1:"))
  vars.filtered %>% filter(Start < BED$Start[1]) %>% select(aachange, Start, End)-> intronic
  print(intronic)
  
  for(i in (seq(1,nrow(BED)-1))){
    print(paste0("Between exon ", i , " and exon " ,i +1,":"))
    vars.filtered %>% filter(Start %in% seq(BED$Stop[i] + 1,BED$Start[i+1] -1)) %>% select(aachange, Start, End)-> intronic
    print(intronic)
  }

  print(paste("After exon", nrow(BED), ":" ))
  vars.filtered %>% filter(Start > BED$Start[nrow(BED)]) %>% select(aachange, Start, End)-> intronic
  print(intronic)


}

test.non.transcript(BED, vars.filtered)


filter.non.transcript<- function(BED, vars.filtered){
  print("Filtering non-transcript changes")
  print(paste("Before exon 1:"))
  vars.filtered %>% filter(Start < BED$Start[1]) %>% select(aachange, Start, End)-> intronic
  vars.filtered %>% filter(!Start < BED$Start[1]) -> res2
  print(intronic)
  
  for(i in (seq(1,nrow(BED)-1))){
    print(paste0("Between exon ", i , " and exon " ,i +1,":"))
    vars.filtered %>% filter(Start %in% seq(BED$Stop[i] + 1,BED$Start[i+1] -1)) %>% select(aachange, Start, End)-> intronic
    res2 %>% filter(!Start %in% seq(BED$Stop[i] + 1,BED$Start[i+1] -1)) -> res2
    print(intronic)
  }
  
  print(paste("After exon", nrow(BED), ":" ))
  vars.filtered %>% filter(Start > BED$Start[nrow(BED)]) %>% select(aachange, Start, End)-> intronic
  res2 %>% filter(!Start > BED$Start[nrow(BED)]) -> res2
  print(intronic)
  return(res2)
  
}

 res <- filter.non.transcript(BED, vars.filtered)




# gene rankings


rank.genes <- function(marv.vars = marv.res) {
#   marv.vars = marv.res
  lof.mutations <- c("frameshift substitution", "frameshift deletion", "frameshift insertion", "stopgain", "stoploss")
  ms.mutations <- c("nonsynonymous SNV", "nonframeshift substitution", "nonframeshift deletion", "nonframeshift insertion")
  marv.vars %>% mutate(ex.funky = sapply(strsplit(levels(marv.res2$ExonicFunc.refGene)[marv.res2$ExonicFunc.refGene],";"), function(x) x[[1]])) -> marv.vars
  marv.vars$mut.category <- "other"
  marv.vars$mut.category[marv.vars$ex.funky %in% lof.mutations] <- "lof"
  marv.vars$mut.category[marv.vars$ex.funky %in% ms.mutations] <- "ms"
  
  ranked.genes <- data.frame(gene = unique(marv.vars$Gene.refGene))
  ranked.genes$marv.lof.count <- NA
  for(i in (1:nrow(ranked.genes))){
    ranked.genes$marv.lof.count[i] <- table(marv.vars$Gene.refGene)[ranked.genes$gene[i]]
  }
#   lof.count <- sapply(ranked.genes$gene, function(x) filter(marv.vars[which(marv.vars$Gene.refGene %in% x),], marv.vars$mut.category == "lof"))
  
#   ranked.genes$lof.count<- lof.count
  
  return(ranked.genes)
  
}
ranked.genes <- rank.genes(marv.res)




head(marv.res)
levels(marv.res$ExonicFunc.refGene)



vec <- c(1,2,3,4)
lapply(vec, function(x) x^2)



a=rep(0:1,5)
b=rep(0,10)
c=rep(1,10)
mat=matrix(cbind(a,b,c),nrow=10,ncol=3)


mat2 <- mat *2 
mat2
