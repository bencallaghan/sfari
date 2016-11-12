# Check Inputs -------------------------------------------------------------

test.fasta <- compare_fastas(translated_cdna = cdna.translated, pp_file = snap2.res, bm_fasta = BM.fasta)
test.files <- check_inputs()
aaDisagreementChecker(vars.filtered,snap2.res) # If aa alignment is low...





# Gene Ranking ------------------------------------------------------------

if(opt.generanks.cache == FALSE){
ranked.genes <- rank_genes2(marv.res3)
ranked.genes <- add_constraint_scores(ranked.genes, exac.constraints)
ranked.list <- actually_rank_genes(ranked.genes)
head(ranked.list, 50)
write.table(ranked.list,"outputs/ranked_list",sep = "\t", row.names = F)           
}

# Metric Correlation ------------------------------------------------------
corr.cadd.snap2 <- cor(vars.filtered$CADD.phred,vars.filtered$PredictProtein,method="spearman")
cor(vars.filtered$cadd,vars.filtered$PredictProtein,method="spearman")

gene.i.scores <- find_gene_scores(vars.filtered, anno.df)

if(!gene.i$name %in% gene.metrics$V1){
write.table(gene.i.scores, path.gene.metrics, append = TRUE)
}
ggplot(sample_n(gene.metrics,1000), aes(x = as.numeric(V1), y = as.numeric(V6), alpha = .5)) + geom_point()


plot_pairs_damage_scores(vars.filtered, dir.cor.plots, as.character(gene.i$name))

#Plot histogram of correlations (genome-wide)
ggplot(gene.metrics.2, aes(x = gs.CADD.v.SNAP2)) + geom_histogram()

# CADD VS SNAP
ggplot(vars.filtered, aes(x = snap2, y = CADD.phred, alpha = .7)) + geom_point()

# Variant Prioritisation --------------------------------------------------
plot_correlation_stuff(vars.filtered,"correlation") #Use these to plot the single correlation plots
plot_correlation_stuff(vars.filtered,"snap2")
plot_correlation_stuff(vars.filtered,"cadd")


prioritised.variants <- prioritise_variants(vars.filtered,lit.variants,marv.res3)
ggplot(prioritised.variants,aes(x=aapos,y=CADD.phred,colour = prioritised, size = !is.na(prioritised)) ) + geom_point()

a<- select_controls_by_aapos2(vars.filtered, c(432, 448))

# SYNGAP1 -----------------------------------------------------------------

head(vars.filtered)

vars.filtered %>% mutate(dist.cs = sqrt((CADD.phred - snap2)^2)) -> vars.filtered
vars.filtered %>% filter(dist.cs > 90) -> high.dist
hist(high.dist$aapos)
ggplot(high.dist, aes(x = snap2, y = CADD.phred)) + geom_point()
vars.filtered %>% filter(dist.cs < 5)  -> low.dist
hist(low.dist$aapos)
ggplot(low.dist, aes(x = snap2, y = CADD.phred)) + geom_point()


cons.res<-  c(432,448,451,465,472,483,484,488,489,525,528,543,562,565,594,595,596,600,601,602,605,615,621,623,625,628,631,634,635,639,650,663,692,693)
cons.res <- c(484, 485, 595,596, 597, 628)
# domain.cons <- select_controls_by_aapos2(vars.filtered,cons.res)
vars.filtered %>% filter(aapos %in% cons.res) %>% dplyr::select(Chr, Start, End, Ref, Alt) -> domain.cons


cadd.quants <- quantile(vars.filtered$CADD.phred, probs = seq(0,1,.1))
snap2.quants <- quantile(vars.filtered$snap2, probs = seq(0,1,.1))

vars.filtered %>% mutate(prioritised = NA) -> pri.vars
marv.vars %>% filter(gene == as.character(gene.i$name)) -> marv.vars

pri.vars %>% filter(exac03 > 0, CADD.phred < cadd.quants[2], snap2 < snap2.quants[2]) -> negs
pri.vars %>% filter(exac03 == 0, CADD.phred > cadd.quants[10], snap2 > snap2.quants[10]) -> pos

pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(negs,2,3,4,5,6)] <- "low-predicted-damage"
pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(pos,2,3,4,5,6)] <- "high-predicted-damage"
pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(lit.vars,1,2,3,4,5)] <- "literature"
pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(marv.vars,3,4,5,6,7)] <- "literature:marv"
pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(domain.cons,1,2,3,4,5)] <- "domain-vars"

##TODO: literature control chooser - just take the position and find the min/max combined damage for each literature variant
#   lit.aas$
# lit.cons <- select_controls_by_aapos(pri.vars, lit.variants$aas)

pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(lit.cons,1,2,3,4,5)] <- "literature.control"

ggplot(pri.vars,aes(x=aapos,y=CADD.phred,colour = prioritised, size = !is.na(prioritised)) ) + geom_point()

##TODO: domain control chooser


#   pri.vars <- pri.vars[which(!is.na(pri.vars$prioritised)),]
return(pri.vars)

