# clean.R:
# Prep / clean all inputs for use
#
# TODO: Incorporate real CADD from kathryn's answers
# Fill in rest of code
# Biomart Transcript ------------------------------------------------------

# BMtranscript %>% filter(transcript_status == "KNOWN", transcript_biotype == "protein_coding",refseq_mrna !="") %>% arrange(desc(transcript_length)) -> BMtranscript.sort
# ensembl_ID <- BMtranscript.sort$ensembl_transcript_id[1]
# refseq_ID <- BMtranscript.sort$refseq_mrna[1]

# Biomart BED -------------------------------------------------------------

BED <- BMbed
cdna_length <- sum(BED$Stop - BED$Start)

# Biomart cDNA ------------------------------------------------------------

CDNA <- BMcdna$cdna
cdna <- tolower(CDNA)
# nchar(CDNA)
cdna.translated <- translate_cdna(cdna)
cdna.orf <- get_cdna_orf(cdna)
orf.coords <- get_orf_coords(cdna)

# Biomart Peptide Sequence ------------------------------------------------

BM.refseqid <- BMpeptide$refseq_peptide
BM.fasta <- BMpeptide$peptide
BM.fasta == cdna.translated
mutfile <- make_mutation_file(BM.fasta)
write.table(mutfile,file = paste0(dir.temp,'mutationfile'))

# Biomart Gene sequence ---------------------------------------------------

genomic_dna <- BMgene$transcript_exon_intron

# Snap2 + All Annovar Info---------------------------------------------------------------

# add cadd.phred
annovar.res$CADD.phred <- sapply(strsplit(as.character(annovar.res2$V2), ","), function(x) x[2]) # add cadd phred scores to annovar results

vars.filtered <- filterGenomicVariants(annovar.res,gene.i$name,gene.i$transcript) #Filter for correct isoform and exonic variants

# Add coordinate string column and cdna positions
vars.filtered$coordinate.string <- coordinate_strings(vars.filtered, 1, 2, 3, 4, 5)
cdnargx <- paste0(".*",gene.i$name,":",gene.i$transcript,":exon[0-9]+:c.([A-Z][0-9]+[A-Z]):p.[A-Z][0-9]+[A-Z].*")
vars.filtered$cdna <- gsub(cdnargx,"\\1",vars.filtered$AAChange.refGene)
vars.filtered$cdna.pos <- gsub("[A-Z]+([0-9]*)[A-Z]+","\\1",vars.filtered$cdna)

# vars.filtered$CADD.phred <- phredScale2(vars.filtered$cadd)
snap2.res$V3 <- gsub("[A-Z]([0-9]+)[A-Z]","\\1",snap2.res$V1)
vars.filtered <- mergePredictProtein(vars.filtered,snap2.res)
vars.filtered <- unfactorizeVariantCols(vars.filtered)
vars.filtered$snap2 <- vars.filtered$PredictProtein

# filter for transcript variants ONLY

vars.filtered <- filter.non.transcript(BED, vars.filtered)


# Gene level Metrics ------------------------------------------------------

gene.metrics %>% filter(!duplicated(gene.metrics$V1)) %>% filter(V2 > -50) %>% arrange(V1) -> gene.metrics.2
colnames(gene.metrics.2) <- c('gene.name', 'gs.gene.length', 'gs.num.aas','gs.num.vars', 'gs.total.vars', 'gs.CADD.v.SNAP2', 
                              'gs.CADD.v.SIFT','gs.CADD.v.PolyPhen2','gs.SNAP2.v.SIFT', 'gs.SNAP2.v.PolyPhen2', 'gs.SIFT.v.PolyPhen2')

# Query variant cleaning  -------------------------------------------------

query.variants <- clean_query_variants(query.variants)

# Input Checks ------------------------------------------------------------



