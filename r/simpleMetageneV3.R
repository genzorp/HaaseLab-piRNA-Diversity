simpleMetagene3 <- function(
                            ## INPUT
                            GR = NULL,
                            RANGE.NAME = NULL,
                            BSSPECIES = NULL,
                            
                            ## OPTIONS
                            USE.READS = TRUE,
                            ALIGN.END = NULL,
                            EXPAND.BY = NULL,
                            
                            ## OUTPUT OPTIONS
                            RETURN.TWO.TABLES = TRUE,
                            
                            ## SETTINGS
                            ALPHABET = c("T","C","G","A")){
  
  ## PG; 09.18.18; Version 1; From small RNA tables V4 
  ## ALIGNS EITHER 5 OR 3 PRIME END
  ## PG; 11.16.19; Version 2; Refined and corrected calculations
  ## PG; 02.03.20; Version 3; BIG CORRECTION
  ## THE ALIGNMENTS WERE CORRECTED
  
  
  ## LIBRARIES
  suppressWarnings(suppressMessages(library("data.table")))
  suppressWarnings(suppressMessages(library("ShortRead")))
  suppressWarnings(suppressMessages(library("GenomicRanges")))
  suppressWarnings(suppressMessages(library("BSgenome")))
  suppressWarnings(suppressMessages(library("hiReadsProcessor")))
  suppressWarnings(suppressMessages(library("BSgenome.Dmelanogaster.UCSC.dm6")))
  suppressWarnings(suppressMessages(library("BSgenome.Mmusculus.UCSC.mm10")))
  suppressWarnings(suppressMessages(library("BSgenome.Hsapiens.UCSC.hg38")))
  
  ## INPUT CHECKING
  if(is.null(GR)){stop("Provide filtered read GR !")}
  if(is.null(RANGE.NAME)){stop("Provide file name to RANGE.NAME !")}
  if(is.null(BSSPECIES)){stop("Provide BSpecies name !!!")}
  if(is.null(EXPAND.BY)){stop("Provide distance value to EXPANDBY !")}
  if(is.null(ALIGN.END)){stop("Provide the end of small RNA to align. '5' or '3' !")}
  
  message("Processing ...")
  message(paste0(" sample: ",RANGE.NAME))
  
  ## EXTEND BOUNDARIES
  expandLimits <- c(-EXPAND.BY, EXPAND.BY)
  
  if(ALIGN.END == 5){
    ## #1 is first nucleotide
    message(" aligning to 5` end")
    GR <- GenomicRanges::resize(GR, width = EXPAND.BY, fix = "start")
    GR <- GenomicRanges::resize(GR, width = EXPAND.BY*2, fix = "end") }
  
  if(ALIGN.END == 3){
    ## #1 is +1 nucleotide
    message(" aligning to 3` end")
    GR <- GenomicRanges::resize(GR, width = EXPAND.BY, fix = "end")
    GR <- GenomicRanges::resize(GR, width = EXPAND.BY*2, fix = "start") }
  
  message(" removing unused contigs")
  GNM <- BSgenome::getBSgenome(eval(parse(text = BSSPECIES)))
  GNM.GR <- GRanges(seqnames = names(GNM), 
                    ranges = IRanges(start = 1, end = GenomeInfoDb::seqlengths(GNM)))
  GR <- IRanges::subsetByOverlaps(GR, GNM.GR, type = "within")
  
  message(" getting sequences (slow)")
  GRS <- BSgenome::getSeq(eval(parse(text = BSSPECIES)), GR)
  
  if(isTRUE(USE.READS)){
    message(" expanding sequences to reads")
    GRS <- replicateReads(GRS, mcols(GR)[["MULT"]]) }
  
  message(" acquiring counts")
  GR.CNT <- as.data.table(ShortRead::alphabetByCycle(GRS, alphabet = ALPHABET))
  GR.CNT[["NAME"]] <- paste("count",RANGE.NAME,gsub("T","U",ALPHABET), sep = "__")
  data.table::setcolorder(GR.CNT, c("NAME",colnames(GR.CNT)[!colnames(GR.CNT) %in% "NAME"]))
  
  message(" calculating frequencies")
  GR.FREQ <- data.table(NAME = GR.CNT[["NAME"]],
                        GR.CNT[,(.SD / lapply(.SD,sum))*100, 
                               .SDcols = colnames(GR.CNT)[!colnames(GR.CNT) %in% "NAME"] ])
  GR.FREQ[["NAME"]] <- gsub("count","frequency",GR.FREQ[["NAME"]])
  
  message(" combining and labeling")
  GR.META <- data.table(data.table::rbindlist(list(GR.CNT,GR.FREQ)))
  colnames(GR.META) <- c("NAME",as.character(c(seq(expandLimits[1],-1),seq(1,expandLimits[2]))))

  
  ## RETURN
  message("Done.")
  
  if(isTRUE(RETURN.TWO.TABLES)){
    TWO.TABLES.L <- list(GR.META[grep("count",GR.META[["NAME"]]),],GR.META[grep("frequency",GR.META[["NAME"]]),])
    names(TWO.TABLES.L) <- c("count","frequency")
    return(TWO.TABLES.L) }
  else{ return(GR.META) }
}