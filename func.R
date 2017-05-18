#####################################
## Functions 
#####################################

####### 1: Gene Finding Utilities -------------------------------------------------------------
findPotentialStartsAndStops <- function(sequence)
{
  # FIND starts and stops
  # sequence <- "aaaaaatga"
  # Define a vector with the sequences of potential start and stop codons
  codons            <- c("atg", "taa", "tag", "tga")
  # Find the number of occurrences of each type of potential start or stop codon
  for (i in 1:4)
  {
    codon <- codons[i]
    # Find all occurrences of codon "codon" in sequence "sequence"
    occurrences <- matchPattern(codon, sequence)
    # Find the start positions of all occurrences of "codon" in sequence "sequence"
    codonpositions <- attr(attr(occurrences,"range"),"start")
    # Find the total number of potential start and stop codons in sequence "sequence"
    numoccurrences <- length(codonpositions)
    if (i == 1)
    {
      # Make a copy of vector "codonpositions" called "positions"
      positions <- codonpositions
      # Make a vector "types" containing "numoccurrences" copies of "codon"
      types <- rep(codon, numoccurrences)
    }
    else
    {
      # Add the vector "codonpositions" to the end of vector "positions":
      positions   <- append(positions, codonpositions, after=length(positions))
      # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
      types       <- append(types, rep(codon, numoccurrences), after=length(types))
    }
  }
  # Sort the vectors "positions" and "types" in order of position along the input sequence:
  indices <- order(positions)
  positions <- positions[indices]
  types <- types[indices]
  # Return a list variable including vectors "positions" and "types":
  mylist <- list(positions,types)
  return(mylist)
}

# 
# findPotentialStartsAndStops <- function(sequence, codons = c("atg", "taa", "tag", "tga")){
#   # find the starting position of the matches for all the codons
#   startpositions <- sapply(codons, function(x){start(matchPattern(x, sequence))})
#   # combine them into a vector and sort by position
#   output <- sort(do.call(c, startpositions))
#   # If we really *need* a list
#   output <- list(as.numeric(output), names(output))
#   return(output)
# }


# Find ORFS
findORFsinSeq <- function(sequence)
{
  require(Biostrings)
  # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
  mylist <- findPotentialStartsAndStops(sequence)
  positions <- mylist[[1]]
  types <- mylist[[2]]
  # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
  orfstarts <- numeric()
  orfstops <- numeric()
  # Make a vector "orflengths" to store the lengths of the ORFs
  orflengths <- numeric()
  # Print out the positions of ORFs in the sequence:
  # Find the length of vector "positions"
  numpositions <- length(positions)
  # There must be at least one start codon and one stop codon to have an ORF.
  if (numpositions >= 2)
  {
    for (i in 1:(numpositions-1))
    {
      posi <- positions[i]
      typei <- types[i]
      found <- 0
      while (found == 0)
      {
        for (j in (i+1):numpositions)
        {
          posj  <- positions[j]
          typej <- types[j]
          posdiff <- posj - posi
          posdiffmod3 <- posdiff %% 3
          # Add in the length of the stop codon
          orflength <- posj - posi + 3
          if (typei == "atg" && (typej == "taa" || typej == "tag" || typej == "tga") && posdiffmod3 == 0)
          {
            # Check if we have already used the stop codon at posj+2 in an ORF
            numorfs <- length(orfstops)
            usedstop <- -1
            if (numorfs > 0)
            {
              for (k in 1:numorfs)
              {
                orfstopk <- orfstops[k]
                if (orfstopk == (posj + 2)) { usedstop <- 1 }
              }
            }
            if (usedstop == -1)
            {
              orfstarts <- append(orfstarts, posi, after=length(orfstarts))
              orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
              orflengths <- append(orflengths, orflength, after=length(orflengths))
            }
            found <- 1
            break
          }
          if (j == numpositions) { found <- 1 }
        }
      }
    }
  }
  # Sort the final ORFs by start position:
  indices <- order(orfstarts)
  orfstarts <- orfstarts[indices]
  orfstops <- orfstops[indices]
  # Find the lengths of the ORFs that we have
  orflengths <- numeric()
  numorfs <- length(orfstarts)
  for (i in 1:numorfs)
  {
    orfstart <- orfstarts[i]
    orfstop <- orfstops[i]
    orflength <- orfstop - orfstart + 1
    orflengths <- append(orflengths,orflength,after=length(orflengths))
  }
  mylist <- list(orfstarts, orfstops, orflengths)
  return(mylist)
}

plotORFsinSeq <- function(sequence)
{
  # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
  mylist <- findPotentialStartsAndStops(sequence)
  positions <- mylist[[1]]
  types <- mylist[[2]]
  # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
  orfstarts <- numeric()
  orfstops <- numeric()
  # Make a vector "orflengths" to store the lengths of the ORFs
  orflengths <- numeric()
  # Print out the positions of ORFs in the sequence:
  numpositions <- length(positions) # Find the length of vector "positions"
  # There must be at least one start codon and one stop codon to have an ORF.
  if (numpositions >= 2)
  {
    for (i in 1:(numpositions-1))
    {
      posi <- positions[i]
      typei <- types[i]
      found <- 0
      while (found == 0)
      {
        for (j in (i+1):numpositions)
        {
          posj <- positions[j]
          typej <- types[j]
          posdiff <- posj - posi
          posdiffmod3 <- posdiff %% 3
          orflength <- posj - posi + 3 # Add in the length of the stop codon
          if (typei == "atg" && (typej == "taa" || typej == "tag" || typej == "tga") && posdiffmod3 == 0)
          {
            # Check if we have already used the stop codon at posj+2 in an ORF
            numorfs <- length(orfstops)
            usedstop <- -1
            if (numorfs > 0)
            {
              for (k in 1:numorfs)
              {
                orfstopk <- orfstops[k]
                if (orfstopk == (posj + 2)) { usedstop <- 1 }
              }
            }
            if (usedstop == -1)
            {
              orfstarts <- append(orfstarts, posi, after=length(orfstarts))
              orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
              orflengths <- append(orflengths, orflength, after=length(orflengths))
            }
            found <- 1
            break
          }
          if (j == numpositions) { found <- 1 }
        }
      }
    }
  }
  # Sort the final ORFs by start position:
  indices <- order(orfstarts)
  orfstarts <- orfstarts[indices]
  orfstops <- orfstops[indices]
  # Make a plot showing the positions of ORFs in the input sequence:
  # Draw a line at y=0 from 1 to the length of the sequence:
  x <- c(1,nchar(sequence))
  y <- c(0,0)
  plot(x, y, ylim=c(0,3), type="l", axes=FALSE, xlab="Nucleotide", ylab="Reading frame", main="Predicted ORFs")
  segments(1,1,nchar(sequence),1)
  segments(1,2,nchar(sequence),2)
  # Add the x-axis at y=0:
  axis(1, pos=0)
  # Add the y-axis labels:
  text(0.9,0.5,"+1")
  text(0.9,1.5,"+2")
  text(0.9,2.5,"+3")
  # Make a plot of the ORFs in the sequence:
  numorfs <- length(orfstarts)
  for (i in 1:numorfs)
  {
    orfstart <- orfstarts[i]
    orfstop <- orfstops[i]
    remainder <- (orfstart-1) %% 3
    if    (remainder == 0) # +1 reading frame
    {
      rect(orfstart,0,orfstop,1,col="cyan",border="black")
    }
    else if (remainder == 1)
    {
      rect(orfstart,1,orfstop,2,col="cyan",border="black")
    }
    else if (remainder == 2)
    {
      rect(orfstart,2,orfstop,3,col="cyan",border="black")
    }
  }
}

translate_cdna <- function(cdnaseq){
  forwardORFs <- findORFsinSeq(cdna)
  plotORFsinSeq(cdna)
  cdna.rev <- c2s(rev(comp(s2c(cdna))))
  reverseORFs <- findORFsinSeq(cdna.rev)
  plotORFsinSeq(cdna.rev)
  if(max(forwardORFs[[3]])>max(reverseORFs[[3]])){
    start_trans <- forwardORFs[[1]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    stop_trans <- forwardORFs[[2]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    cdna.translated <-  c2s(seqinr::translate(s2c(myorf)))                           
  }else{
    stop_trans <- reverseORFs[[1]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    stop_trans <- reverseORFs[[2]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    cdna.translated <- c2s(seqinr::translate(s2c(myorf)))
    
  }
  return(cdna.translated)
}

get_cdna_orf <- function(cdnaseq){
  forwardORFs <- findORFsinSeq(cdna)
  plotORFsinSeq(cdna)
  cdna.rev <- c2s(rev(comp(s2c(cdna))))
  reverseORFs <- findORFsinSeq(cdna.rev)
  plotORFsinSeq(cdna.rev)
  if(max(forwardORFs[[3]])>max(reverseORFs[[3]])){
    start_trans <- forwardORFs[[1]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    stop_trans <- forwardORFs[[2]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    #     cdna.translated <-  c2s(seqinr::translate(s2c(myorf)))                           
  }else{
    stop_trans <- reverseORFs[[1]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    stop_trans <- reverseORFs[[2]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    #     cdna.translated <- c2s(seqinr::translate(s2c(myorf)))
    
  }
  return(myorf)
}

genomic_dna_for_annovar <- function(genomic_dna, bed){
  print(genomic_dna)
  #   str(genomic_dna)
  print(paste0("length: ",nchar(genomic_dna)))
  print(bed)
  
  if(bed$Strand[1] == 1){
    gStart <- bed$Start[1]
    gStop <- bed$Stop[nrow(bed)]
    Chrom <- bed$Chrom[1]
  }else{#Testing for minus strand genes here...
    gStart <- bed$Start[nrow(bed)]
    gStop <- bed$Stop[1]
    Chrom <- bed$Chrom[1]
  }
  
  if(nchar(genomic_dna) == gStop + 1 - gStart){
    print("Genomic DNA length agrees with BED coordinates, calculating annovar file....")
    annovar.df <- data.frame(Chrom = Chrom,
                             Start = gStart:gStop,
                             Stop = gStart:gStop,
                             Ref = s2c(genomic_dna))
    
    annovar.df <- annovar.df[rep(1:nrow(annovar.df),each=4),]
    annovar.df$Alt <- rep(c("A","C","G","T"),nchar(genomic_dna))
    annovar.df %>% filter(Alt != Ref) -> annovar.df
    return(annovar.df)
  } else{
    print("Genomic DNA length does not agree with BED coordinates... fix before continuing")
  }
}
# anno.df <- genomic_dna_for_annovar(BMgene$transcript_exon_intron,BED)
# head(anno.df,20)

get_orf_coords <- function(cdnaseq){
  forwardORFs <- findORFsinSeq(cdna)
  plotORFsinSeq(cdna)
  cdna.rev <- c2s(rev(comp(s2c(cdna))))
  reverseORFs <- findORFsinSeq(cdna.rev)
  plotORFsinSeq(cdna.rev)
  if(max(forwardORFs[[3]])>max(reverseORFs[[3]])){
    start_trans <- forwardORFs[[1]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    stop_trans <- forwardORFs[[2]][which(forwardORFs[[3]] == max(forwardORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    #     cdna.translated <-  c2s(seqinr::translate(s2c(myorf)))                           
  }else{
    stop_trans <- reverseORFs[[1]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    stop_trans <- reverseORFs[[2]][which(reverseORFs[[3]] == max(reverseORFs[[3]]))]
    myorf <- substring(cdna, start_trans, stop_trans)
    #     cdna.translated <- c2s(seqinr::translate(s2c(myorf)))
    
  }
  return(c(start_trans,stop_trans))
}

###### 2:  Functions for Dealing with other programs-----------------------------------------------------------

scribe_awkcommand <- function(CHROM,BED,GENE,INPUTDIR){
  awkcommand <- paste0("awk '$1 == ",CHROM, " && (")
  for (i in 1:nrow(BED)){
    if(i < nrow(BED)){
      print(i)
      exoniccmd <- (paste0("($2 > ",BED$Start[i] - 1, " && $2 < ", BED$Stop[i] + 1 ,") || "))
      awkcommand <- paste0(awkcommand,exoniccmd)
    }else{
      print(i)
      exoniccmd <- (paste0("($2 > ",BED$Start[i] - 1, " && $2 < ", BED$Stop[i] + 1,")) "))
      awkcommand <- paste0(awkcommand,exoniccmd,"{print $0}' /space/bin/annovar/humandb/hg19_cadd.txt > ",
                           INPUTDIR,"cadd",GENE, " 2> ", INPUTDIR, "cadd", GENE, ".err" )
    }
  }
  return(awkcommand)
}

make_anno_file <- function(cdna.orf){
  for(i in 1:length(s2c(cdna.orf))){
    print(s2c(cdna.orf)[i])
  }
  anno.df <- data.frame(matrix(nrow=nchar(cdna.orf),ncol=5))
  colnames(anno.df) <- c("CHR","START","STOP","REF","ALT")
  anno.df$REF <- s2c(cdna.orf)
}

make_mutation_file <- function(FASTA){
  FASTA <- s2c(FASTA)
  aas <- AA_STANDARD
  FASTA <- FASTA[which(FASTA %in% aas)]
  plength <- length(FASTA)
  mutfile <- vector('character')
  for(i in 1:plength){
    for(j in 1:20){
      mut <- paste0(FASTA[i],i,aas[j])
      mutfile<- append(mutfile,mut)
    }
  }
  return(mutfile)
}

compare_fastas <- function(translated_cdna, pp_file, bm_fasta){
  x <- gsub("([A-Z][0-9]+)[A-Z]","\\1",levels(pp_file$V1)[pp_file$V1],perl=TRUE) 
  x <- x[which(duplicated(x) == FALSE)]
  x <- paste(gsub("([A-Z])[0-9]+","\\1",x, perl = TRUE),collapse="", sep = "")
  pp_fasta <- paste0(x,"*")
  
  print(paste0("Biomart FASTA, refseq ID ", BM.refseqid))
  print(bm_fasta)
  
  print(paste0("Translated cdna FASTA, refseq ID ", refseq_ID))
  print(translated_cdna)
  
  print("PredictProtein Parsed FASTA")
  print(pp_fasta)
  
  if(pp_fasta == translated_cdna && translated_cdna == bm_fasta){
    print("Input and calculated fastas all agree")
    return(TRUE)
  } else {
    print("*** Input and calculated fastas do not agree, fix before continuing ***")
    return(FALSE)
  }
}

phredScale <- function(scores){
  scores.scaled <- scale(scores,)
  scores.prob <- dnorm(scores.scaled)
  phred <- c(-10*log10(scores.prob))
  phred[which(phred > 40)] <- 40
  return(phred)
  
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

phredScale2 <- function(scores){
  phred <- -10*log10(1-rank(scores)/max(rank(scores)))
  phred[!is.finite(phred)] <- 40
  return(phred)
}
# phredScale2(1:1000)
# phredScale2(gene_vars_merge$cadd)



check_inputs <- function(){
  miss <- NULL
  if(file.exists(path.fasta) == FALSE){
    miss <- "fasta"
  } 
  #   if(file.exists(bedpath) == FALSE){
  #     miss <- c(miss,"bed")
  #   } 
  #   if(file.exists(caddpath) == FALSE){
  #     miss <- c(miss,"cadd")
  
  
  #   } 
  if(file.exists(paste0(dir.inputs,gene.i$name,"_anno.hg19_multianno.csv")) == FALSE){
    miss <- c(miss,"caddanno")
  } 
  if(is.null(miss)){
    print("**All inputs exist, continuing with analysis**")
    return(TRUE)
  }else{
    print(paste("**",miss, "files missing, fix before prioritisation**"))
    return(FALSE)
  }
}

filterGenomicVariants <- function(genomic_variants,gene,refseq_transcript){
  head(genomic_variants)
  genomic_variants$Func.refGene <- gsub("exonic;splicing","exonic",genomic_variants$Func.refGene)
  genomic_variants %>% filter(grepl("exonic", Func.refGene)) -> filtered_variants
  filtered_variants %>% filter(grepl(paste0(gene,":",refseq_transcript,".+"),AAChange.refGene)) -> filtered_variants
  
  aachange_regex <- paste0(".*",gene,":",refseq_transcript,":exon[0-9]+:c.[A-Z][0-9]+[A-Z]:p.([A-Z][0-9]+[A-Z]).*")
  filtered_variants$aachange <- gsub(aachange_regex,"\\1",filtered_variants$AAChange.refGene)
  
  aapos_regex <- "[A-Z]([0-9]+)[A-Z]"
  filtered_variants$aapos <- as.numeric(gsub(aapos_regex,"\\1",filtered_variants$aachange))
  return(filtered_variants)
}

mergePredictProtein <- function(genomic_variants,predict_protein){
  cat(nrow(predict_protein)/19, "Predict Protein amino acids\n")
  cat(nrow(genomic_variants)/9, "Isoform amino acids\n")
  cat(sum(genomic_variants$aachange %in% predict_protein$V1)/nrow(genomic_variants) * 100,"% Amino acid Alignment\n")# Predict Protein doesn't score stopgains/losses
  
  colnames(predict_protein) <- c("aachange.pp","PredictProtein","aapos.pp")
  variants_merged <- merge(genomic_variants, predict_protein, by.x="aachange", by.y="aachange.pp",all.x=TRUE)
  
  variants_merged %>% arrange((Start)) -> variants_merged
  return(variants_merged)
}

aaDisagreementChecker <- function(genomic_variants,predict_protein){
  cat("ALL==============")
  cat(genomic_variants$aachange[which(!genomic_variants$aachange %in% predict_protein$V1)])
  cat(length(genomic_variants$aachange[which(!genomic_variants$aachange %in% predict_protein$V1)]))
  cat("SYNONYMOUS===============")
  #   x <- gsub(([A-Z])[0-9]+([A-Z]),)
}

unfactorizeVariantCols <- function(genomic_variants){
  genomic_variants$CADD.phred <- as.numeric(genomic_variants$CADD.phred)
  genomic_variants$cadd <- as.numeric(genomic_variants$cadd)
  genomic_variants$PredictProtein <- as.numeric(genomic_variants$PredictProtein)
  genomic_variants$exac03 <- as.numeric(as.character(genomic_variants$exac03))
  genomic_variants$exac03[is.na(genomic_variants$exac03)] <- 0
  return(genomic_variants)
}

df <- data.frame('num' <- c(1,2,3),"alpha" <- c('a','b','c'))
unfactorizeColumns <- function(df,columns){
  df <- df
  columns
  for(i in columns){
    print(i)
    
  }
}


bioGetTranscript<- function(GENE){
  attributes <- c("refseq_mrna","ensembl_transcript_id","transcript_start","transcript_end" ,
                  "transcript_status","transcript_count", "transcript_biotype", "transcript_source",
                  "transcript_version","transcript_length")
  BMtranscript <- getBM(attributes = attributes, 
                        filters = c("hgnc_symbol"),values = list(GENE), mart = mart, verbose = FALSE)
  BMtranscript %>% filter(transcript_status == "KNOWN", transcript_biotype == "protein_coding",refseq_mrna !="") %>% 
    arrange(desc(transcript_length)) -> BMtranscript.sort
}
bioGetBED <- function(GENE){
  
}

#_---3-3-3-3-3-3-3-3--3-3-3-3-3-3-3-3-3-3--3-3-3
retrieve_exonic_vcf <- function(GENE){
  # Take in gene, match to Biomart
  # Find canonical transcript (transcript of largest size)
  BMtranscript.sort <- bioGetTranscript(GENE)
  ensembl_ID <- BMtranscript.sort$ensembl_transcript_id[1]
  refseq_ID <- BMtranscript.sort$refseq_mrna[1]
  if(nrow(BMtranscript.sort) == 0){
    return(FALSE) # If not in BM, exit early and return blank
  }
  # Find Chromosome
  # Retrieve exonic BED for canonical transcript
  BED <- getBM(attributes = c("chromosome_name","exon_chrom_start","exon_chrom_end" ,"rank","strand","ensembl_transcript_id"), 
               filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(GENE,ensembl_ID), mart = mart, verbose = FALSE)
  CHROM <- BED$chromosome_name[1]
  colnames(BED) <- c("Chrom", "Start", "Stop", "rank", "strand", "ensembl_transcript_id")
  BED$Gene <- GENE
  cdna_length <- sum(BED$Stop - BED$Start)
  # Get exonic DNA sequence
  BMgene <- getBM(attributes = c("transcript_exon_intron","refseq_mrna","ensembl_transcript_id"), 
                  filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(GENE,ensembl_ID), mart = mart, verbose = FALSE)
  BMgenecoords <- getBM(attributes = c("start_position","end_position","genomic_coding_start","genomic_coding_end", "ensembl_transcript_id"), 
                        filters = c("hgnc_symbol","ensembl_transcript_id"),values = list(GENE,ensembl_ID), mart = mart, verbose = FALSE)
  
  #   cdna.orf
  genomic_dna <- BMgene$transcript_exon_intron
  # Build VCF
  vcf <- genomic_dna_for_annovar(BMgene$transcript_exon_intron,BED)
  # Get CDNA sequence
  return(vcf)
}


merge_multianno_and_snap2 <- function(multiannodf, snap2df){
  # Merges annovar output and snap2 results
  # Filters for exonic variants in canonical transcript (by filterGenomicVariants())
  # Sanity checks for length and mapping of amino acids
  print("Merging annovar and snap2 results...")
  head(multiannodf)
  head(snap2df)
  
  multiannodf <- filterGenomicVariants(multiannodf,GENE,TRANSCRIPT)
  snap2df$V1 <- levels(snap2df$V1)[snap2df$V1]
  
  cat(nrow(snap2df)/19, "Predict Protein amino acids\n")
  cat(snap2df$V1[1], "...",snap2df$V1[nrow(snap2df)],"\n")
  cat(nrow(multiannodf)/9, "Isoform amino acids\n")
  cat(multiannodf$aachange[1], "...",multiannodf$aachange[nrow(multiannodf)],"\n")
  cat(sum(multiannodf$aachange %in% snap2df$V1)/nrow(multiannodf)*100,"% Amino acid Alignment\n")# Predict Protein doesn't score stopgains/losses
  
  multiannodf$CADD.phred <- phredScale2(multiannodf$cadd)
  snap2df$V3 <- gsub("[A-Z]([0-9]+)[A-Z]","\\1",snap2df$V1)
  
  colnames(snap2df) <- c("snap2aachange","snap2","snap2aapos")
  mergedDF <- merge(multiannodf, snap2df, by.x="aachange", by.y="snap2aachange") 
  mergedDF %>% arrange((Start)) -> mergedDF
  
  return(mergedDF)
}



plot_correlation_stuff <- function(variants, mode = "correlation"){
  corr <- function(variants){
    outpath <- paste0(dir.outputs, as.character(gene.i$name), "CADD_SNAP2_Correlation.png")
    p1 <- ggplot(variants,aes(x= CADD.phred, y = snap2)) + geom_point()
    ggsave(outpath,p1)
  }
  snap2 <- function(variants){
    outpath <- paste0(dir.outputs, as.character(gene.i$name), "SNAP2_Scores_byAAPOS.png")
    p1 <- ggplot(variants,aes(x= aapos, y = snap2)) + geom_point()
    ggsave(outpath,p1)
  }
  cadd <- function(variants){
    outpath <- paste0(dir.outputs, as.character(gene.i$name), "CADD_SCORES_byAAPOS.png")
    p1 <- ggplot(variants,aes(x= aapos, y = CADD.phred)) + geom_point()
    ggsave(outpath, p1)
  }
  switch(mode, 
         correlation = corr(variants),
         snap2 = snap2(variants),
         cadd = cadd(variants)
  )
  
}

plot_pairs_damage_scores <- function(multianno.df, outfolder, genename){
  multianno.df %>% dplyr::select(snap2, CADD.phred, LJB_PolyPhen2, LJB_SIFT, exac03, aapos) -> damage.scores
  damage.scores$PolyPhen2 <- as.numeric(levels(damage.scores$LJB_PolyPhen2)[damage.scores$LJB_PolyPhen2])
  damage.scores$SIFT <- as.numeric(levels(damage.scores$LJB_SIFT)[damage.scores$LJB_SIFT])
  damage.scores$Snap2 <- as.numeric(damage.scores$snap2)
#   damage.scores$exac03 <- as.numeric(levels(damage.scores$exac03)[damage.scores$exac03])
  damage.scores <- mutate(damage.scores, in.exac = ifelse(!is.na(exac03), TRUE, FALSE))
  damage.scores$in.exac <- as.factor(damage.scores$in.exac)
  damage.scores %>% dplyr::select(aapos, Snap2, SIFT, PolyPhen2, CADD.phred, in.exac) -> damage.scores
  
  pairsplot <- ggpairs(damage.scores, aes(colour = in.exac, alpha = 0.4),upper = list(
    continuous = wrap('cor', method = "spearman")))
  
  #   ggpairs(damage.scores, aes(colour = in.exac, alpha = 0.4),upper = list(
  #     continuous = wrap('cor', method = "spearman")),lower = list(
  #       continuous = wrap('cor', method = "pearson")))
  
  # png(file.path(outfolder, paste0(genename, ".png")), height=1000, width=1000)
  # print(pairsplot)
  # dev.off()
  outpath <- paste0(outfolder, as.character(gene.i$name), "_gene_metrics_pairsplot.png")
  ggsave(outpath, pairsplot)
  
  
}
# plot_pairs_damage_scores(multianno.2, "/home/bcallaghan/NateDBCopy/20kcorr/pairsplots/", GENE)


find_gene_scores <- function(multianno.df, vcf){
  # make annotated df work nicely, return all the gene level info we have
  multianno.df %>% dplyr::select(snap2, CADD.phred, LJB_PolyPhen2, LJB_SIFT, exac03, aapos) -> damage.scores
  damage.scores$PolyPhen2 <- as.numeric(levels(damage.scores$LJB_PolyPhen2)[damage.scores$LJB_PolyPhen2])
  damage.scores$SIFT <- as.numeric(levels(damage.scores$LJB_SIFT)[damage.scores$LJB_SIFT])
  damage.scores$Snap2 <- as.numeric(damage.scores$snap2)
#   damage.scores$exac03 <- as.numeric(levels(damage.scores$exac03)[damage.scores$exac03])
  damage.scores <- mutate(damage.scores, in.exac = ifelse(!is.na(exac03), TRUE, FALSE))
  damage.scores$in.exac <- as.factor(damage.scores$in.exac)
  damage.scores %>% dplyr::select(aapos, Snap2, SIFT, PolyPhen2, CADD.phred, in.exac) -> damage.scores
  
  gs.gene.length <- abs(max(vcf$Start) - min(vcf$Stop))
  gs.num.aas <- max(damage.scores$aapos)
  gs.num.vars <- nrow(damage.scores)
  gs.total.vars <- nrow(vcf)
  gs.CADD.v.SNAP2 <- cor(damage.scores$CADD.phred, damage.scores$Snap2, method="spearman",use="pairwise.complete.obs")
  gs.CADD.v.SIFT <- cor(damage.scores$CADD.phred, damage.scores$SIFT, method="spearman",use="pairwise.complete.obs")
  gs.CADD.v.PolyPhen2 <- cor(damage.scores$CADD.phred, damage.scores$PolyPhen2, method="spearman",use="pairwise.complete.obs")
  
  gs.SNAP2.v.SIFT <- cor(damage.scores$Snap2, damage.scores$SIFT, method="spearman",use="pairwise.complete.obs")
  gs.SNAP2.v.PolyPhen2 <- cor(damage.scores$Snap2, damage.scores$PolyPhen2, method="spearman",use="pairwise.complete.obs")
  
  gs.SIFT.v.PolyPhen2 <- cor(damage.scores$SIFT, damage.scores$PolyPhen2, method="spearman",use="pairwise.complete.obs")
  
  gene.info.string <- paste(as.character(gene.i$name), gs.gene.length, gs.num.aas,gs.num.vars, gs.total.vars, gs.CADD.v.SNAP2, 
                            gs.CADD.v.SIFT,gs.CADD.v.PolyPhen2,gs.SNAP2.v.SIFT, gs.SNAP2.v.PolyPhen2, 
                            gs.SIFT.v.PolyPhen2, sep="\t")
  gene.info.vec <- c(as.character(gene.i$name), gs.gene.length, gs.num.aas,gs.num.vars, gs.total.vars, gs.CADD.v.SNAP2, 
                            gs.CADD.v.SIFT,gs.CADD.v.PolyPhen2,gs.SNAP2.v.SIFT, gs.SNAP2.v.PolyPhen2, 
                            gs.SIFT.v.PolyPhen2)
  return(gene.info.vec)
}

check_empty_gene_error <- function(vcf){
  if(is.logical(vcf) == TRUE){
    return(TRUE)
  }else{
    FALSE
  }
}

#3 - Gene level functions - ranking, correlations, etc ----------------------

rank.genes <- function(marv.vars = marv.res) {
  
  lof.mutations <- c("frameshift substitution", "frameshift deletion", "frameshift insertion", "stopgain", "stoploss")
  ms.mutations <- c("nonsynonymous SNV", "nonframeshift substitution", "nonframeshift deletion", "nonframeshift insertion")
  marv.vars %>% mutate(ex.funky = sapply(strsplit(levels(marv.vars$ExonicFunc.refGene)[marv.vars$ExonicFunc.refGene],";"), function(x) x[[1]])) -> marv.vars
  marv.vars$mut.category <- "other"
  marv.vars$mut.category[marv.vars$ex.funky %in% lof.mutations] <- "lof"
  marv.vars$mut.category[marv.vars$ex.funky %in% ms.mutations] <- "ms"
  ranked.genes <- data.frame(gene = unique(marv.vars$Gene.refGene))
  ranked.genes$marv.count <- sapply(ranked.genes$gene, function(x) nrow(marv.vars[which(marv.vars$Gene.refGene %in% x),]))
  ranked.genes$marv.lof.count <- sapply(ranked.genes$gene, function(x) nrow(marv.vars[which(marv.vars$Gene.refGene %in% x),] %>% filter(mut.category == 'lof')))
  ranked.genes$marv.ms.count <- sapply(ranked.genes$gene, function(x) nrow(marv.vars[which(marv.vars$Gene.refGene %in% x),] %>% filter(mut.category == 'ms')))
  ranked.genes <- ranked.genes[order(ranked.genes$marv.count,decreasing=TRUE),]
  return(ranked.genes)
}

count_marv_variants <- function(marv.vars, method = "denovo"){
  lof.mutations <- c("frameshift substitution", "frameshift deletion", "frameshift insertion", "stopgain", "stoploss")
  ms.mutations <- c("nonsynonymous SNV", "nonframeshift substitution", "nonframeshift deletion", "nonframeshift insertion")
  marv.vars$mut.category <- "other"
  marv.vars$mut.category[marv.vars$category %in% lof.mutations] <- "lof"
  marv.vars$mut.category[marv.vars$category %in% ms.mutations] <- "ms"
  if(method == "denovo"){
    marv.vars %>% filter(denovo == "yes") -> marv.vars
  }
  if(method == "total"){
    marv.vars %>% filter(denovo == "yes" | denovo == "unknown" | denovo == "NULL") -> marv.vars
  }
  ranked.genes <- data.frame(gene = unique(marv.vars$gene))
  ranked.genes$marv.count <- sapply(ranked.genes$gene, function(x) nrow(marv.vars[which(marv.vars$gene %in% x),]))
  ranked.genes$marv.lof.count <- sapply(ranked.genes$gene, function(x) nrow(marv.vars[which(marv.vars$gene %in% x),] %>% filter(mut.category == 'lof')))
  ranked.genes$marv.ms.count <- sapply(ranked.genes$gene, function(x) nrow(marv.vars[which(marv.vars$gene %in% x),] %>% filter(mut.category == 'ms')))
  ranked.genes <- ranked.genes[order(ranked.genes$marv.count,decreasing=TRUE),]
  return(ranked.genes)
}

add_constraint_scores <- function(ranked.genes, exac.constraints){
  ranked.genes$MSZ <- exac.constraints[match(ranked.genes$gene, exac.constraints$gene, nomatch=NA), c("mis_z")]
  ranked.genes$pLI <- exac.constraints[match(ranked.genes$gene, exac.constraints$gene, nomatch=NA), c("pLI")]
  ranked.genes$pRec <- exac.constraints[match(ranked.genes$gene, exac.constraints$gene, nomatch=NA), c("pRec")]
  ranked.genes$pNull <- exac.constraints[match(ranked.genes$gene, exac.constraints$gene, nomatch=NA), c("pNull")]
  ranked.genes
}

range01 <- function(x,...){(x-min(x,...))/(max(x,...)-min(x,...))}


actually_rank_genes <- function(ranked.genes){
  # ranked.genes <- ranked.genes[!is.na(ranked.genes$MSZ),]
  ranked.genes$MSZ[is.na(ranked.genes$MSZ)] <- 0
  ranked.genes$pLI[is.na(ranked.genes$pLI)] <- 0
  
  ranked.genes %>% 
  mutate(MSZ = replace(MSZ, MSZ <= 0, 0)) %>%
  mutate(scale.msz = (MSZ-min(MSZ))/(max(MSZ)-min(MSZ))) %>%
  mutate(lof.ms.aggregate = range01(marv.ms.count * scale.msz + marv.lof.count * pLI, na.rm = T) * 100 ) %>% 
  arrange(desc(lof.ms.aggregate)) %>% 
  dplyr::select(gene,marv.count, marv.lof.count, marv.ms.count, MSZ, pLI, pRec, pNull, scale.msz, lof.ms.aggregate) -> ranked.list
  return(ranked.list)
  
}


gene_level_correlations <- function(df){
  df %>% filter(!is.na(gs.gene.length)) %>% 
    dplyr::select(gs.gene.length,gs.num.aas,gs.num.vars,gs.CADD.v.SNAP2,gs.CADD.v.SIFT,
                  gs.CADD.v.PolyPhen2, gs.SNAP2.v.SIFT, gs.SNAP2.v.PolyPhen2, gs.SIFT.v.PolyPhen2) -> df
  
  head(df)
  ggpairs(df,aes(alpha = .1))
}
# gene_level_correlations(x)

# messy_genes <- function(df, features){
#   for(i in features){
#     print(features[i])
#     df %>% dplyr::select(i) %>% quantile(na.rm = TRUE, 0.05)
#   }
#   quantile(x$gs.CADD.v.SNAP2, na.rm = TRUE, .05)
#   
# }

messy_genes <- function(df, percentile){
  featurenames <- c('gs.CADD.v.SNAP2','gs.CADD.v.SIFT','gs.CADD.v.PolyPhen2', 'gs.SNAP2.v.SIFT', 'gs.SNAP2.v.PolyPhen2', 'gs.SIFT.v.PolyPhen2')
  thresh <- quantile(df$gs.CADD.v.SNAP2, na.rm = TRUE, percentile)
  print(thresh)
  df %>% filter(gs.CADD.v.SNAP2 < thresh) %>% dplyr::select(GENE) -> a
  
  thresh <- quantile(df$gs.CADD.v.SIFT, na.rm = TRUE, percentile)
  print(thresh)
  df %>% filter(gs.CADD.v.SIFT < thresh) %>% dplyr::select(GENE) -> b
  
  thresh <- quantile(df$gs.CADD.v.PolyPhen2, na.rm = TRUE, percentile)
  print(thresh)
  df %>% filter(gs.CADD.v.PolyPhen2 < thresh) %>% dplyr::select(GENE) -> c
  
  thresh <- quantile(df$gs.SNAP2.v.SIFT, na.rm = TRUE, percentile)
  print(thresh)
  df %>% filter(gs.SNAP2.v.SIFT < thresh) %>% dplyr::select(GENE) -> d
  
  thresh <- quantile(df$gs.SNAP2.v.PolyPhen2, na.rm = TRUE, percentile)
  print(thresh)
  df %>% filter(gs.SNAP2.v.PolyPhen2 < thresh) %>% dplyr::select(GENE) -> e
  
  thresh <- quantile(df$gs.SIFT.v.PolyPhen2, na.rm = TRUE, percentile)
  print(thresh)
  df %>% filter(gs.SIFT.v.PolyPhen2 < thresh) %>% dplyr::select(GENE) -> f
  
  max.len <- max(length(a),length(b),length(c),length(d),length(e),length(f))
  #   a = c(a, rep(NA, max.len - length(a)))
  #   b = c(b, rep(NA, max.len - length(b)))
  #   c = c(c, rep(NA, max.len - length(c)))
  #   d = c(d, rep(NA, max.len - length(d)))
  #   e = c(e, rep(NA, max.len - length(e)))
  #   f = c(f, rep(NA, max.len - length(f)))
  
  #   res <- data.frame(c(a,b,c,d,e,f))
  print(paste("Gene names < ", percentile, "correlation:", featurenames[1]))
  print(a)
  print(paste("Gene names < ", percentile, "correlation:", featurenames[2]))
  print(b)
  print(paste("Gene names < ", percentile, "correlation:", featurenames[3]))
  print(c)
  print(paste("Gene names < ", percentile, "correlation:", featurenames[4]))
  print(d)
  print(paste("Gene names < ", percentile, "correlation:", featurenames[5]))
  print(e)
  print(paste("Gene names < ", percentile, "correlation:", featurenames[6]))
  print(f)
  
  a <- as.character(a$GENE)
  b <- as.character(b$GENE)
  c <- as.character(c$GENE)
  d <- as.character(d$GENE)
  e <- as.character(e$GENE)
  f <- as.character(f$GENE)
  
  res <- Reduce(intersect, list(a,b,c,d,e,f))
  return(res)
  #   return(res)
  #   quantile(x$gs.CADD.v.SNAP2, na.rm = TRUE, .05)
}
# genelist <- messy_genes(x, 0.05)

filter.non.transcript<- function(BED, vars.filtered){
  print("Filtering non-transcript changes")
  print(paste("Before exon 1:"))
  vars.filtered %>% filter(Start < BED$Start[1]) %>% dplyr::select(aachange, Start, End)-> intronic
  vars.filtered %>% filter(!Start < BED$Start[1]) -> res2
  print(intronic)
  
  for(i in (seq(1,nrow(BED)-1))){
    print(paste0("Between exon ", i , " and exon " ,i +1,":"))
    vars.filtered %>% filter(Start %in% seq(BED$Stop[i] + 1,BED$Start[i+1] -1)) %>% dplyr::select(aachange, Start, End)-> intronic
    res2 %>% filter(!Start %in% seq(BED$Stop[i] + 1,BED$Start[i+1] -1)) -> res2
    print(intronic)
  }
  
  print(paste("After exon", nrow(BED), ":" ))
  vars.filtered %>% filter(Start > BED$Start[nrow(BED)]) %>% dplyr::select(aachange, Start, End)-> intronic
  res2 %>% filter(!Start > BED$Stop[nrow(BED)]) -> res2
  print(intronic)
  return(res2)
  
}

prioritise_variants <- function(variants,lit.vars,marv.vars){
  cadd.quants <- quantile(variants$CADD.phred, probs = seq(0,1,.1))
  snap2.quants <- quantile(variants$CADD.phred, probs = seq(0,1,.1))
  
  vars.filtered %>% mutate(prioritised = NA) -> pri.vars
  marv.vars %>% filter(gene == as.character(gene.i$name)) -> marv.vars
  
  pri.vars %>% filter(exac03 > 0, CADD.phred < cadd.quants[2], snap2 < snap2.quants[2]) -> negs
  pri.vars %>% filter(exac03 == 0, CADD.phred > cadd.quants[10], snap2 > snap2.quants[10]) -> pos
  
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(negs,2,3,4,5,6)] <- "low-predicted-damage"
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(pos,2,3,4,5,6)] <- "high-predicted-damage"
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(lit.vars,1,2,3,4,5)] <- "literature"
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(marv.vars,3,4,5,6,7)] <- "literature:marv"
  
  
  ##TODO: literature control chooser - just take the position and find the min/max combined damage for each literature variant
  pri.vars %>% filter(prioritised %in% c("literature", "literature:marv")) %>% dplyr::select(aapos) -> lit.aas
  lit.cons <- select_controls_by_aapos(pri.vars, lit.aas)
  
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(lit.cons,1,2,3,4,5)] <- "literature.control"
  
  
  ##TODO: domain control chooser
  
  
#   pri.vars <- pri.vars[which(!is.na(pri.vars$prioritised)),]
  return(pri.vars)
}

prioritise_variants2 <- function(variants,lit.vars,marv.vars){
  cadd.quants <- quantile(variants$CADD.phred, probs = seq(0,1,.1))
  snap2.quants <- quantile(variants$CADD.phred, probs = seq(0,1,.1))
  
  vars.filtered %>% mutate(prioritised = NA) -> pri.vars
  marv.vars %>% filter(gene == as.character(gene.i$name)) -> marv.vars
  
  pri.vars %>% filter(exac03 > 0, CADD.phred < cadd.quants[2], snap2 < snap2.quants[2]) -> negs
  pri.vars %>% filter(exac03 == 0, CADD.phred > cadd.quants[10], snap2 > snap2.quants[10]) -> pos
  
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(negs,2,3,4,5,6)] <- "low-predicted-damage"
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(pos,2,3,4,5,6)] <- "high-predicted-damage"
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(lit.vars,1,2,3,4,5)] <- "literature"
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(marv.vars,3,4,5,6,7)] <- "literature:marv"
  
  
  ##TODO: literature control chooser - just take the position and find the min/max combined damage for each literature variant
#   lit.aas$
  lit.cons <- select_controls_by_aapos(pri.vars, lit.aas)
  
  pri.vars$prioritised[coordinate_strings(pri.vars,2,3,4,5,6) %in% coordinate_strings(lit.cons,1,2,3,4,5)] <- "literature.control"
  
  
  ##TODO: domain control chooser
  
  
  #   pri.vars <- pri.vars[which(!is.na(pri.vars$prioritised)),]
  return(pri.vars)
}



select_controls_by_aapos <- function(pri.vars, aapos){
  df <- data.frame(Chr= numeric(0), Start= numeric(0), End = numeric(0), Ref = character(0), Alt = character(0))
  for(i in 1:length(lit.aas$aapos)){
#     print(i)
    pri.vars %>% filter(aapos == lit.aas$aapos[i]) %>% filter(!prioritised %in% c('literature', 'literature:marv')) %>% filter(CADD.phred == max(CADD.phred)) -> max.damage
    pri.vars %>% filter(aapos == lit.aas$aapos[i]) %>% filter(!prioritised %in% c('literature', 'literature:marv')) %>% filter(CADD.phred == min(CADD.phred)) -> min.damage
    controls <- dplyr::select(rbind(max.damage,min.damage),c(Chr, Start, End, Ref, Alt))
    df <- rbind(df,controls)       
  }
  return(df)
}

select_controls_by_aapos2 <- function(pri.vars, aapos){
  df <- data.frame(Chr= numeric(0), Start= numeric(0), End = numeric(0), Ref = character(0), Alt = character(0))
  for(i in 1:length(aapos)){
    #     print(i)
    pri.vars %>% filter(aapos == aapos[i])%>% filter(CADD.phred == max(CADD.phred)) -> max.damage
    pri.vars %>% filter(aapos == aapos[i]) %>% filter(CADD.phred == min(CADD.phred)) -> min.damage
    controls <- dplyr::select(rbind(max.damage,min.damage),c(Chr, Start, End, Ref, Alt))
    df <- rbind(df,controls)       
  }
  return(df)
}

select_controls_by_domain <- function(pri.vars, aapos){
  df <- data.frame(Chr= numeric(0), Start= numeric(0), End = numeric(0), Ref = character(0), Alt = character(0))
  for(i in 1:length(lit.aas$aapos)){
    #     print(i)
    pri.vars %>% filter(aapos == lit.aas$aapos[i]) %>% filter(!prioritised %in% c('literature', 'literature:marv')) %>% filter(CADD.phred == max(CADD.phred)) -> max.damage
    pri.vars %>% filter(aapos == lit.aas$aapos[i]) %>% filter(!prioritised %in% c('literature', 'literature:marv')) %>% filter(CADD.phred == min(CADD.phred)) -> min.damage
    controls <- dplyr::select(rbind(max.damage,min.damage),c(Chr, Start, End, Ref, Alt))
    df <- rbind(df,controls)       
  }
  return(df)
}
# a <- select_controls_by_aapos(prioritised.variants,lit.aas)


coordinate_strings <- function(df,chrnum,startnum, stopnum,refnum,altnum){
  #chrnum, startnum etc correspond to the column numbers containing said information
  # Returns a vector of strings of length nrow(df) 
  coord.strng <- paste0(df[,chrnum],":",df[,startnum],"-",df[,stopnum],":",df[,refnum],'/' ,df[,altnum])
  #coordinate_strings(vars.filtered,2,3,4,5,6)
  return(coord.strng)
}

search_variants <- function(variants, query.variants, mode = "coordinates"){
  if(mode == coordinatestrin){
    
  }
}

# Query functions --------------------------------
# Query all gene variants with subset list 

clean_query_variants <- function(query.variants){#*make this a funtion
  for(i in 1:nrow(query.variants)){
    print(i)
    if(query.variants$querytype[i] == "coordinates"){
    query.variants$coordinate.string[i] <- coordinate_strings(query.variants[i,], 1,2,2,3,4)
    }else{
    query.variants$coordinate.string[i] <- NA
    }
  }
  return(query.variants)
}

query_variants <- function(variants, query){
  #   res <- data.frame()
  res <- data.frame(0, matrix(nrow = nrow(query), ncol = ncol(vars.filtered)))
  colnames(res) <- colnames(vars.filtered)
  for(i in 1:nrow(res)){ #nrow(res)
    print(paste0("query", i))
    if(query$querytype[i] == "coordinates"){
      print("Query type: coordinates")
      if(query$coordinate.string[i] %in% variants$coordinate.string){
        res[i,] <- lapply(variants[match(as.character(query$coordinate.string[i]) ,as.character(variants$coordinate.string)),], as.character)
      } else {
        res$coordinate.string[i] <- query$coordinate.string[i]
      }
      res$Source[i] <- as.character(query$Source[i])
      res$denovo[i] <- as.character(query$denovo[i])
    }
    
    if(query$querytype[i] == "cdna"){
      print("Query type: cdna")
      if(query$cdna[i] %in% variants$cdna){
        res[i,] <- lapply(variants[match(as.character(query$cdna[i]) ,as.character(variants$cdna)),],as.character)
      } else {
        res$cdna[i] <- as.character(query$cdna[i])
      }
      res$Source[i] <- as.character(query$Source[i])
      res$denovo[i] <- as.character(query$denovo[i])      
      
    }
    if(query$querytype[i] == "protein"){
      print("Query type: aachange")
      if(query$aachange[i] %in% variants$aachange){
        res[i,] <- lapply(variants[match(as.character(query$aachange[i]) ,as.character(variants$aachange)),], as.character)
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

logfile_writeline <- function(line){
  cat(paste0(line,'\n'), file = log_con)
}

