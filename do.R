# Check Inputs -------------------------------------------------------------

test.fasta <- compare_fastas(translated_cdna = cdna.translated, pp_file = snap2.res, bm_fasta = BM.fasta)
test.files <- check_inputs()
aaDisagreementChecker(vars.filtered,snap2.res) # If aa alignment is low...





# Gene Ranking ------------------------------------------------------------

# First attempt
marv.res

marv.res2 <- marv.res
marv.res2$ex.funky <- sapply(strsplit(levels(marv.res2$ExonicFunc.refGene)[marv.res2$ExonicFunc.refGene],";"), function(x) x[[1]])
# gene.ranking <- data.frame(gene = unique(marv.res2$Gene.refGene), lof.counts)

ranked.genes <- rank.genes(marv.res)

ranked.genes$MSZ <- exac.constraints[match(ranked.genes$gene, exac.constraints$gene, nomatch=NA), c("mis_z")]
ranked.genes$pLI <- exac.constraints[match(ranked.genes$gene, exac.constraints$gene, nomatch=NA), c("pLI")]
ranked.genes$pRec <- exac.constraints[match(ranked.genes$gene, exac.constraints$gene, nomatch=NA), c("pRec")]
ranked.genes$pNull <- exac.constraints[match(ranked.genes$gene, exac.constraints$gene, nomatch=NA), c("pNull")]

ranked.genes %>% 
#   mutate(lof.score = pLI * marv.lof.count) %>% 
  filter(MSZ > -10) %>%
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
write.table(ranked.list,"outputs/ranked_list",sep = "\t")           
ranked.genes %>% arrange(desc(lof.ms.score)) %>% head(20)
ranked.genes %>% arrange(desc(lof.score)) %>% head(20)
ranked.genes %>% arrange(desc(ms.score)) %>% head(20)
# write.table(x=ranked.genes,"outputs/ranked_genes_constraint",sep="\t", row.names=F)
head(marv.res2)

#Second attempt
marv.res3




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



# SYNGAP1 -----------------------------------------------------------------

head(vars.filtered)

vars.filtered %>% mutate(dist.cs = sqrt((CADD.phred - snap2)^2)) -> vars.filtered
vars.filtered %>% filter(dist.cs > 90) -> high.dist
hist(high.dist$aapos)
ggplot(high.dist, aes(x = snap2, y = CADD.phred)) + geom_point()
vars.filtered %>% filter(dist.cs < 5)  -> low.dist
hist(low.dist$aapos)
ggplot(low.dist, aes(x = snap2, y = CADD.phred)) + geom_point()







