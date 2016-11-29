
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



ranked.genes <- rank_genes2(marv.res3)
ranked.genes <- add_constraint_scores(ranked.genes, exac.constraints)


ranked.genes %>% 
  #   mutate(lof.score = pLI * marv.lof.count) %>% 
  filter(MSZ > -10) %>%
  replace(MSZ <= 0 , 0) %>%
  #   mutate(ms.score = MSZ * marv.ms.count) %>% 
  #   mutate(lof.ms.score = rank(ms.score * lof.score)) %>%
  #   mutate(rank.ms = rank(marv.ms.count)) %>%
  #   mutate(rank.lof = rank(marv.lof.count)) %>%
  mutate(scale.msz = (MSZ-min(MSZ))/(max(MSZ)-min(MSZ))) %>%
  mutate(lof.ms.aggregate = (marv.ms.count * scale.msz + marv.lof.count * pLI)) %>% 
  #   arrange(desc(marv.count)) %>%
  #   arrange(desc(lof.ms.score)) %>% head(20)
  arrange(desc(lof.ms.aggregate)) %>% 
  select(gene,marv.count, marv.lof.count, marv.ms.count, MSZ, pLI, pRec, pNull, scale.msz, lof.ms.aggregate) -> ranked.list
ranked.list %>% head(100)




# 04-10-16 SFARI variants -------------------------------------------------


a <- read.table("/home/bcallaghan/Projects/Onenineteen/119_dn_variants2")

# Create an annovar input file
write.table(a, file = paste0("/home/bcallaghan/Projects/Onenineteen/119varsyeahyeah"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)

# Run annovar a second time for Phred-scaled CADD scores
if(opt.annovar.cache & file.exists(path.annovar.out2)){
  print("Using cached annovar output2")
}else{
  annocmd2 <- paste0("perl /space/bin/annovar/annotate_variation.pl ", "/home/bcallaghan/Projects/Onenineteen/119varsyeahyeah ", " /space/bin/annovar/humandb -filter -dbtype cadd -buildver hg19 -out ",
                     "/home/bcallaghan/Projects/Onenineteen/119varsyeahyeahres", "_anno2 -otherinfo")
  cmd.out <- run.remote(cmd=annocmd2 , remote= "-q apu",stderr.redirect=F)
}


# SYNGAPVARIANTS ----------------------------------------------------------

# query.variants <- read.table("inputs/SYNGAP_query_variants.csv", fill=T,header = T, sep = "\t")
# query.variants <- read.table("inputs/SYNGAP_query_variants2.csv", fill=T,header = T, sep = ",")
# vars.filtered$cdna.pos <- gsub("[A-Z]+([0-9]*)[A-Z]+","\\1",vars.filtered$cdna)

# for(i in 1:nrow(query.variants)){
#   print(i)
#   if(query.variants$querytype[i] == "coordinates"){
#     query.variants$coordinate.string[i] <- coordinate_strings(query.variants[i,], 1,2,2,3,4)
#   }else{
#     query.variants$coordinate.string[i] <- NA
#   }
# }
vars.filtered2 <- as.character(vars.filtered)

variants <- vars.filtered #Testing
query <- query.variants #testing
# i = 65
query_variants <- function(variants, query){
#   res <- data.frame()
  res <- data.frame(0, matrix(nrow = nrow(query), ncol = ncol(vars.filtered)))
  colnames(res) <- colnames(vars.filtered)
  for(i in 1:nrow(res)){ #nrow(res)
    print(paste0("query", i))
    if(query$querytype[i] == "coordinates"){
      print("Query type coords")
      if(query$coordinate.string[i] %in% variants$coordinate.string){
        res[i,] <- variants[match(as.character(query$coordinate.string[i]) ,as.character(variants$coordinate.string)),]
      } else {
        res$coordinate.string[i] <- query$coordinate.string[i]
      }
      res$Source[i] <- as.character(query$Source[i])
      res$denovo[i] <- as.character(query$denovo[i])
    }
    
    if(query$querytype[i] == "cdna"){
      print("query type cdna")
      if(query$cdna[i] %in% variants$cdna){
        res[i,] <- variants[match(as.character(query$cdna[i]) ,as.character(variants$cdna)),]
      } else {
        res$cdna[i] <- as.character(query$cdna[i])
      }
      res$Source[i] <- as.character(query$Source[i])
      res$denovo[i] <- as.character(query$denovo[i])      
      
    }
    if(query$querytype[i] == "protein"){
      print("query type aachange")
      if(query$aachange[i] %in% variants$aachange){
        res[i,] <- variants[match(as.character(query$aachange[i]) ,as.character(variants$aachange)),]
      } else {
        res$aachange[i] <- as.character(query$aachange[i])
      }
      res$Source[i] <- as.character(query$Source[i])
      res$denovo[i] <- as.character(query$denovo[i])      
      
    }
    if(query$querytype[i] == "site"){
      print("site")
      sitevars <- variants[which(variants$aapos %in% query$aachange[i]),]
      maxsnap <- sitevars[which(sitevars$snap2 %in% max(sitevars$snap2))[1],]
      res[i,] <- maxsnap
      res$Source[i] <- as.character(query$Source[i])
      res$denovo[i] <- as.character(query$denovo[i])    
    }
  }
  return(res)
}

queried <- query_variants(vars.filtered, query.variants)
queried

queried %>% select(aachange, cdna, coordinate.string, Func.refGene, CADD.phred, snap2, exac03,Source, denovo ) -> res.queried

vars.filtered %>% 
  filter(snap2 < 0, CADD.phred < 10, exac03 > 0) %>% 
  select(aachange, cdna, coordinate.string, Func.refGene, CADD.phred, snap2, exac03) %>% 
  mutate(Source = "negative controls") %>% mutate(denovo = NA) -> res.negs

8%>% 
  select(aachange, cdna, coordinate.string, Func.refGene, CADD.phred, snap2, exac03) %>% 
  mutate(Source = "positive controls") %>% mutate(denovo = NA) -> res.pos ; res.pos


# cdna.sub.start <- 490
# cdna.sub.stop <- 490
# substr(cdna, cdna.sub.start, cdna.sub.stop )
# marv.res3 %>% filter(gene == "SYNGAP1")
write.table(x=rbind(res.queried,res.negs), file=paste0("outputs/SYNGAP1_variants_",format(Sys.time(), '%m_%d_%H.%M'),".csv"), sep = ",", row.names=FALSE)

# filtered.list <- query_by_coordinates(vars.filtered, query.variants)



