miRbase2BED <- function(miRBASEFILE = NULL,
                        saveBed = FALSE,
                        miRNAOnly = TRUE){
  
  ## Pavol Genzor; 
  ## 08.08.18; Version 1
  ## 12.18.19; Version 2; Update
  
  ## LIBRARIES
  suppressWarnings(suppressMessages(library("data.table")))
  
  ## INPUT CHECKING
  if(is.null(miRBASEFILE)){stop("Please provide .gtf2 or .gff3 file from miRbase")}
  
  ## DETECT FILE TYPE
  FILE.TYPE <- tstrsplit(tail(unlist(tstrsplit(miRBASEFILE, split = "/")), n=1), split = "\\.")[[2]]
  
  if(FILE.TYPE %in% "gff3"){
    message("Processing .gff3 file")
    ## READ
    GFF <- fread(miRBASEFILE, header = FALSE)

    ## ADD COLUMN NAMES AND SELECT USEFUL ONES
    GFF.COLNAMES <- c("chr","score0","feature","start","end","score1","strand","score2","description")
    colnames(GFF) <- GFF.COLNAMES

    ## SPLIT DESCRIPTION AND TAKE OU THE NAME AND THE ORIGIN
    ID = tstrsplit(tstrsplit(GFF$description, split = ";")[[1]], split = "=")[[2]]
    ALIAS = tstrsplit(tstrsplit(GFF$description, split = ";")[[2]], split = "=")[[2]]
    NAME = tstrsplit(tstrsplit(GFF$description, split = ";")[[3]], split = "=")[[2]]
    ORIGIN = tstrsplit(tstrsplit(GFF$description, split = ";")[[4]], split = "=")[[2]]
    
    ## MAKE BED TABLE
    miRBED <- data.table(chr = GFF$chr, 
                         start = GFF$start,
                         end = GFF$end,
                         name = NAME,
                         feature = GFF$feature,
                         strand = GFF$strand,
                         gene_id = ID,
                         alias = ALIAS,
                         origin =  ORIGIN)
    
    ## ONLY MATURE MIRNAS
    if(isTRUE(miRNAOnly)){miRBED <- miRBED[miRBED$feature %in% "miRNA",]} }
  
  if(FILE.TYPE %in% "gtf2"){
    message("Processing .gtf2 file")
    ## READ
    GTF2 <- fread(miRBASEFILE, header = FALSE)
    
    ## ADD COLUMN NAMES AND SELECT USEFUL ONES
    GTF2.COLNAMES <- c("chr","database","feature","start","end","score1","strand","score2","description")
    colnames(GTF2) <- GTF2.COLNAMES
    
    ## SPLIT THE DESCRIPTION AND TAKE OUT GENE AND TRANSCRIPT IDs
    GID <- tstrsplit(tstrsplit(tstrsplit(GTF2$description, split = ";")[[1]], split =" ")[[2]], split ="\"")[[2]]
    TID <- tstrsplit(tstrsplit(tstrsplit(GTF2$description, split = ";")[[2]], split =" ")[[3]], split ="\"")[[2]]
    
    ## MAKE BED TABLE
    miRBED <- data.table(chr = GTF2$chr, 
                         start = GTF2$start,
                         end = GTF2$end,
                         gene_id = GID,
                         transript_id = TID,
                         strand = GTF2$strand,
                         database = GTF2$database) }
  
  if(isTRUE(saveBed)){
    message("\tsaving .bed file")
    FILE.NAME <- tstrsplit(tail(unlist(tstrsplit(miRBASEFILE, split = "/")), n=1), split = "\\.")[[1]]
    write.table(x = miRBED, file = paste(FILE.NAME,".bed", sep = ""), 
                quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    message(paste("\tsaved",paste(FILE.NAME,".bed", sep = ""),"file in the working directory", sep = " ")) }
  
  message("Done.")
  return(miRBED)
}