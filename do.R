# Check Inputs -------------------------------------------------------------

test.fasta <- compare_fastas(translated_cdna = cdna.translated, pp_file = snap2.res, bm_fasta = BM.fasta)
test.files <- check_inputs()
aaDisagreementChecker(vars.filtered,snap2.res) # If aa alignment is low 


# Gene Ranking ------------------------------------------------------------

if(opt.generanks.cache == FALSE){
  ranked.genes <- count_marv_variants(marv.res, method = "total")
  # ranked.genes <- count_marv_variants(marv.res, method = "denovo")
  ranked.genes <- add_constraint_scores(ranked.genes, exac.constraints)
  ranked.list <- actually_rank_genes(ranked.genes)
  ranked.list$brain.expressed <- ranked.list$gene %in% brain.expressed.gene.list$V1
  head(ranked.list, 50)
  write.table(ranked.list,"outputs/ranked_list.csv",sep = ",", row.names = F) # write the most up to date in outputs dir
  # write.table(ranked.list, paste0(dir.outputs, "ranked_list.csv"), sep = ",", row.names = F) # Also write in current gene output dir
  
  ranked.genes <- count_marv_variants(marv.res, method = "denovo")
  ranked.genes <- add_constraint_scores(ranked.genes, exac.constraints)
  ranked.list <- actually_rank_genes(ranked.genes)
  ranked.list$brain.expressed <- ranked.list$gene %in% brain.expressed.gene.list$V1
  head(ranked.list, 50)
  write.table(ranked.list,"outputs/ranked_list_denovo.csv",sep = ",", row.names = F) # write the most uptodate in outputs dir
  
  
}

# Metric Correlation ------------------------------------------------------
corr.caddp.snap2 <- cor(as.numeric(vars.filtered$CADD.phred),as.numeric(vars.filtered$snap2,method="spearman"), use = "complete.obs")
cat(corr.caddp.snap2)
corr.cadd.snap2 <- cor(as.numeric(vars.filtered$cadd),as.numeric(vars.filtered$snap2),method="spearman", use = "complete.obs")

gene.i.scores <- find_gene_scores(vars.filtered, anno.df)

if(!gene.i$name %in% gene.metrics$V1){
write.table(gene.i.scores, path.gene.metrics, append = TRUE)
}
ggplot(sample_n(gene.metrics,1000), aes(x = as.numeric(V1), y = as.numeric(V6), alpha = .5)) + geom_point()

# plot_pairs_damage_scores(vars.filtered, dir.cor.plots, as.character(gene.i$name))
# p1.pairs <- plot_pairs_damage_scores(vars.filtered, dir.cor.plots, as.character(gene.i$name))
# p1.pairs

#Plot histogram of correlations (genome-wide)
# ggplot(gene.metrics.2, aes(x = gs.CADD.v.SNAP2)) + geom_histogram()

# CADD VS SNAP
# ggplot(vars.filtered, aes(x = snap2, y = CADD.phred, alpha = .7)) + geom_point()

# Variant Prioritisation --------------------------------------------------
##### Replaced by query vvvv, slated for removal
# prioritised.variants <- prioritise_variants(vars.filtered,lit.variants,marv.res)
# ggplot(prioritised.variants,aes(x=aapos,y=CADD.phred,colour = prioritised, size = !is.na(prioritised)) ) + geom_point()


# Literature / biochemical assay variant query ---------------------------
print("Querying for variants of interest:")
print(query.variants)
variants <- vars.filtered; query <- query.variants
queried <- query_variants(vars.filtered, query.variants)
queried %>% select(aachange, cdna, coordinate.string, Func.refGene, ExonicFunc.refGene, CADD.phred, snap2, exac03,Source, denovo ) -> res.queried
res.queried


# Calculated control variants --------------------------------------------
cadd.quants <- quantile(vars.filtered$CADD.phred, probs = seq(0,1,.1), na.rm = T)
snap2.quants <- quantile(vars.filtered$snap2, probs = seq(0,1,.1), na.rm = T)

vars.filtered %>% 
  filter(snap2 < 0, CADD.phred < 15, exac03 > 0) %>% 
  select(aachange, cdna, coordinate.string, Func.refGene, ExonicFunc.refGene, CADD.phred, snap2, exac03) %>% 
  mutate(Source = "negative controls") %>% mutate(denovo = NA) -> res.population.controls


vars.filtered %>% 
  filter(snap2 > 0, CADD.phred > 20, exac03 == 0) %>% 
  select(aachange, cdna, coordinate.string, Func.refGene, ExonicFunc.refGene, CADD.phred, snap2, exac03) %>% 
  mutate(Source = "calculated positive controls") %>% mutate(denovo = NA) -> res.calculated.positives

# Outputs ----------------------------------------
# Output Variant of interest file
# write.table(x=rbind(res.queried,res.negs), file=paste0("outputs/SYNGAP1_variants_",format(Sys.time(), '%m_%d_%H.%M'),".csv"), sep = ",", row.names=FALSE)
write.table(x=rbind(res.queried,res.population.controls), file=paste0(dir.outputs,gene.i$name,"_VariantsOfInterest",".csv"), sep = ",", row.names=FALSE)
write.table(x=rbind(res.calculated.positives), file=paste0(dir.outputs,gene.i$name,"_CaculatedHighImpact",".csv"), sep = ",", row.names=FALSE)

# Output Plots ------------------------------------------------------------
# Plot CADD vs SNAP2 correlation, CADD by AApos, SNAP2 by AApos
# Saves in dir.outputs path
plot_correlation_stuff(vars.filtered,"correlation") #Use these to plot the single correlation plots
plot_correlation_stuff(vars.filtered,"snap2")
plot_correlation_stuff(vars.filtered,"cadd")

#Plot histogram of correlations (genome-wide)
# Saves in dir.outputs path
ggplot(gene.metrics.2, aes(x = gs.CADD.v.SNAP2)) + geom_histogram()
ggsave(paste0(dir.outputs, "20k_snap2_vs_CADD_correlations.png"))

# Plot gene metric pairsplot
# Saves in dir.outputs path
plot_pairs_damage_scores(vars.filtered, dir.outputs, as.character(gene.i$name))




# SYNGAP1 -----------------------------------------------------------------
# 
# head(vars.filtered)
# 
# vars.filtered %>% mutate(dist.cs = sqrt((CADD.phred - snap2)^2)) -> vars.filtered
# vars.filtered %>% filter(dist.cs > 90) -> high.dist
# hist(high.dist$aapos)
# ggplot(high.dist, aes(x = snap2, y = CADD.phred)) + geom_point()
# vars.filtered %>% filter(dist.cs < 5)  -> low.dist
# hist(low.dist$aapos)
# ggplot(low.dist, aes(x = snap2, y = CADD.phred)) + geom_point()
# 
# 
# cons.res<-  c(432,448,451,465,472,483,484,488,489,525,528,543,562,565,594,595,596,600,601,602,605,615,621,623,625,628,631,634,635,639,650,663,692,693)
# cons.res <- c(484, 485, 595,596, 597, 628)
# # domain.cons <- select_controls_by_aapos2(vars.filtered,cons.res)
# vars.filtered %>% filter(aapos %in% cons.res) %>% dplyr::select(Chr, Start, End, Ref, Alt) -> domain.cons
# 
# 
# 
# vars.filtered %>% mutate(prioritised = NA) -> pri.vars
# marv.vars %>% filter(gene == as.character(gene.i$name)) -> marv.vars
# 
# pri.vars %>% filter(exac03 > 0, CADD.phred < cadd.quants[2], snap2 < snap2.quants[2]) -> negs
# pri.vars %>% filter(exac03 == 0, CADD.phred > cadd.quants[10], snap2 > snap2.quants[10]) -> pos
# 
# pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(negs,2,3,4,5,6)] <- "low-predicted-damage"
# pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(pos,2,3,4,5,6)] <- "high-predicted-damage"
# pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(lit.vars,1,2,3,4,5)] <- "literature"
# pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(marv.vars,3,4,5,6,7)] <- "literature:marv"
# pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(domain.cons,1,2,3,4,5)] <- "domain-vars"
# 
# ##TODO: literature control chooser - just take the position and find the min/max combined damage for each literature variant
# #   lit.aas$
# # lit.cons <- select_controls_by_aapos(pri.vars, lit.variants$aas)
# 
# pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(lit.cons,1,2,3,4,5)] <- "literature.control"
# 
# ggplot(pri.vars,aes(x=aapos,y=CADD.phred,colour = prioritised, size = !is.na(prioritised)) ) + geom_point()
# 
# ##TODO: domain control chooser
# 
# 
# #   pri.vars <- pri.vars[which(!is.na(pri.vars$prioritised)),]
# return(pri.vars)

