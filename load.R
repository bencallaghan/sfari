########################################
# Load in all files:
# MARVdb
# BioMart Files: Transcripts, BED, Peptide seq, gene sequence
# Annovar output
# Snap2 
# Query file (for querying literature variants)
#
########################################
 

# Biomart -----------------------------------------------------------------

# Biomart Setup

listMarts(host="grch37.ensembl.org", path="/biomart/martservice" )
mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host="grch37.ensembl.org",path="/biomart/martservice")
BMattr <- listAttributes(mart)
BMattr[grep("position",BMattr$name),]
attributes <- c("refseq_mrna","ensembl_transcript_id","transcript_start","transcript_end" ,
                "transcript_status","transcript_count", "transcript_biotype", "transcript_source",
                "transcript_version","transcript_length")

# Transcript
BMtranscript <- getBM(attributes = attributes, filters = c("hgnc_symbol"),values = list(as.character(gene.i$name)), mart = mart, verbose = FALSE)
# Canonical transcript
BMtranscript %>% filter(transcript_status == "KNOWN", transcript_biotype == "protein_coding",refseq_mrna !="") %>% arrange(desc(transcript_length)) -> BMtranscript.sort

ensembl_ID <- BMtranscript.sort$ensembl_transcript_id[1]
refseq_ID <- BMtranscript.sort$refseq_mrna[1]
# ensembl_ID <- BMtranscript.sort$ensembl_transcript_id[3]
# refseq_ID <- BMtranscript.sort$refseq_mrna[3]
# BED
BMbed <- getBM(attributes = c("chromosome_name","exon_chrom_start","exon_chrom_end" ,"rank","strand","ensembl_transcript_id"), 
               filters = c("chromosome_name","hgnc_symbol","ensembl_transcript_id"),
               values = list(as.character(gene.i$chromosome), as.character(gene.i$name), ensembl_ID), mart = mart, verbose = FALSE)
colnames(BMbed) <- c("Chrom","Start","Stop","Exon","Strand","Transcript")
BMbed$Gene <- as.character(gene.i$name)

# Transcript (cDNA)
BMcdna <- getBM(attributes = c("cdna","refseq_mrna","ensembl_transcript_id"), 
                filters = c("chromosome_name","hgnc_symbol","ensembl_transcript_id"),
                values = list(as.character(gene.i$chromosome), as.character(gene.i$name), ensembl_ID), mart = mart, verbose = FALSE)


# Peptide Sequence --------------------------------------------------------
BMpeptide <- getBM(attributes = c("peptide","refseq_mrna","refseq_peptide","ensembl_transcript_id"), 
                   filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(as.character(gene.i$name),ensembl_ID), mart = mart, verbose = FALSE)

# Genomic Sequence and Coordinates ---------------------------------------------
BMgene <- getBM(attributes = c("transcript_exon_intron","refseq_mrna","ensembl_transcript_id"), 
                filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(as.character(gene.i$name),ensembl_ID), mart = mart, verbose = FALSE)
BMgenecoords <- getBM(attributes = c("start_position","end_position","genomic_coding_start","genomic_coding_end", "ensembl_transcript_id"), 
                      filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(as.character(gene.i$name),ensembl_ID), mart = mart, verbose = FALSE)

# Snap2 -------------------------------------------------------------------

# Snap2 mapping
snap2map <- read.table(paste0(dir.inputs,"fastamap")) # Maps Gene Symbols to snap2 results files
snapfile <- snap2map[which(snap2map$V2 == as.character(as.character(gene.i$name))),1]
snap2path <- paste0('/misc/pipeline42/ppdatabases/snap2results/UP000005640_9606/', snapfile, '.snap2.parsed')


# Snap2 results -----------------------------------------------------------
# If snap2 results aren't coming in (and you're running locally), may need to change server (ex from otto to apu)
# Otherwise, check that pipeline24 is mounted on the rstudio server (and that you can see the snap2 files)

# if(opt.session.local == TRUE){ #At least as long as it's on the rtest... this might work better
#  snap2.res <- read.table(pipe(paste0("ssh -p 22000 echu113@otto.pavlab.chibi.ubc.ca cat ", snap2path))) 
# }else {
snap2.res  <- read.table(snap2path)
# }

# Annovar  -----------------------------------------------------------

anno.df <- genomic_dna_for_annovar(BMgene$transcript_exon_intron,bed = BMbed)

# Create an annovar input file
write.table(anno.df, file = paste0(dir.temp,as.character(gene.i$name),"_anno_in"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)

# Create and run annovar command - build annovar annotation file for gene.i
if(opt.annovar.cache & file.exists(path.annovar.out)){
  print("using cached annovar output")
}else{
annocmd <- paste0("perl /space/bin/annovar/table_annovar.pl ", path.anno.in, " /space/bin/annovar/humandb/ -buildver hg19 -out ", dir.temp,as.character(gene.i$name), "_anno -remove -protocol refGene,genomicSuperDups,esp6500si_all,1000g2012apr_all,snp135,ljb_all,exac03,cadd -otherinfo -operation g,r,f,f,f,f,f,f -nastring . -csvout")
cmd.out <- run.remote(cmd=annocmd , remote= "apu")
}

# Run annovar a second time for Phred-scaled CADD scores
if(opt.annovar.cache & file.exists(path.annovar.out2)){
  print("Using cached annovar output2")
}else{
annocmd2 <- paste0("perl /space/bin/annovar/annotate_variation.pl ", path.anno.in, " /space/bin/annovar/humandb -filter -dbtype cadd -buildver hg19 -out ",
                   dir.temp,as.character(gene.i$name), "_anno2 -otherinfo")
cmd.out <- run.remote(cmd=annocmd2 , remote= "-q apu",stderr.redirect=F)
}
# Annovar Results ---------------------------------------------------------

annovar.res <- read.csv(path.annovar.out) #change to annovar.res
annovar.res2 <- read.table(path.annovar.out2, sep="\t")

# Fasta -------------------------------------------------------------------
# fasta <- read.table(path.fasta)


# Gene Level Metrics ------------------------------------------------------

gene.metrics <- read.table(path.gene.metrics, header = F, sep = "\t", fill=TRUE)

# MARVdb Variants ---------------------------------------------------------

# marv.res <- read.csv(path.marv) # Change to marv.res
# marv.res <- read.table("inputs/marv.allvariants_27-09-2016.csv",sep="\t", header = T)
#marv.res <- read.table("inputs/marvdbdump-03-07-17.csv", header = T, fill = T, sep = "\t") 
#marv.res <- read.table("inputs/marvdbdump-04-04-17.csv", header = T, fill = T, sep = "\t", comment.char = "" ) 
marv.res <- read.table("inputs/marv_variants.csv", header = T, fill = T, sep = "\t", comment.char = "")


# Pull in ranked gene list (from prior runs - check opt.generanks.cache to redo this ranking)

if(opt.generanks.cache == TRUE){
gene.list <- read.table("outputs/ranked_list",sep = "\t", header = T)           
}
# Exac constraints file

exac.constraints <- read.table("inputs/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt", header = T)


# Brain expressed gene list from Marjan
# Should do my own sometime with GTEx data
brain.expressed.gene.list <- read.table("inputs/marjan_brain_expressed_gene_list.txt")


#* Literature variants (in vcf format) --- slated to remove ------
# 
# if(file.exists(path.gene.lit.variants)){
#   lit.variants <- read.table(path.gene.lit.variants,sep = "\t", col.names=c("chr","start","stop","ref","alt","source"),fill=TRUE)
# }

# Query variants

# query.variants <- read.table("inputs/DYRK1A_query_variants.csv", fill=T,header = T, sep = ",") 
query.variants <- read.table(path.query.variants, fill=T,header = T, sep = "," ) 
# query.variants <- read.table("inputs/DYRK1A_query_variants_0413.csv", fill=T,header = T, sep = "," ) 

###### Custom sequence input
# Uncomment if you need to manually add sequences
#
#
#




