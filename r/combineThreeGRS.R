combineThreeGRS <- function(
                            ## INPUT
                            GRL=NULL,
                            REPLICATE.NAMES=NULL,
                            
                            ## SETTINGS
                            MC.CORES=3
                            ){
  
  ## AUTHOR: Pavol Genzor
  ## 11.27.19: Version 1
  
  ## NOTE: This function takes time
  
  ## LIBRARIES
  suppressMessages(suppressWarnings(library("data.table")))
  suppressMessages(suppressWarnings(library("plyr")))
  suppressMessages(suppressWarnings(library("dplyr")))
  suppressMessages(suppressWarnings(library("GenomicRanges")))
  suppressMessages(suppressWarnings(library("stringi")))
  suppressMessages(suppressWarnings(library("parallel")))
  
  ## INPUT CHECKING
  if(is.null(GRL)) stop("Provide named list of Genomic Ranges")
  if(isFALSE(class(GRL) %in% c("list"))) stop("Input must be a list of GRanges")
  if(isTRUE(is.null(names(GRL)))) stop("This list needs to be named")
  if(is.null(REPLICATE.NAMES)) stop("Provide names for three replicates")
  if(!length(REPLICATE.NAMES)==3) stop("Provide exactly THREE replicates")
  
  message("Proceeding ...")
  
  message(" combining replicates")
  T.REPS.GR <- c(GRL[[REPLICATE.NAMES[1]]],
                 GRL[[REPLICATE.NAMES[2]]],
                 GRL[[REPLICATE.NAMES[3]]])
  
  message(" making a data.table")
  T.REP.DT <- as.data.table(T.REPS.GR)
  
  message(" removing width")
  T.REP.DT[[c("width")]] <- NULL
  
  message(" creating a unique ID")
  T.REP.DT[["id"]] <- stringi::stri_c(T.REP.DT[["seqnames"]],
                                      T.REP.DT[["start"]], 
                                      T.REP.DT[["end"]], 
                                      T.REP.DT[["strand"]], sep = "__")
  
  message(" summarizing mcols")
  USQ.COUNT <- T.REP.DT[,.N, by = "id"]
  USQ.NH <- T.REP.DT[,lapply(.SD, max), by = "id", .SDcols = c("NH")]
  USQ.MULT <- T.REP.DT[,lapply(.SD, sum), by = "id", .SDcols = c("MULT")]
  
  message(" combining summary results") ## SLOWER
  NF <- setDT(plyr::join_all(dfs = list(USQ.COUNT,USQ.NH,USQ.MULT), 
                             by = "id", type = "full"))
  
  message(" expanding the ID")
  EXPANDED.ID.L <- mclapply(seq(1,4), mc.cores = MC.CORES, function(i){
    dplyr::nth(data.table::tstrsplit(NF[["id"]],split="__",fixed=TRUE),i) })
  names(EXPANDED.ID.L) <- c("chr","start","end","strand")
  
  message(" adding expanded ID items to the table")
  NF[["chr"]] <- EXPANDED.ID.L[["chr"]]
  NF[["start"]] <- EXPANDED.ID.L[["start"]]
  NF[["end"]] <- EXPANDED.ID.L[["end"]]
  NF[["strand"]] <- EXPANDED.ID.L[["strand"]]
  
  #NF[, `:=` (chr=dplyr::nth(data.table::tstrsplit(id,split="__",fixed=TRUE),1),
  #           start=dplyr::nth(data.table::tstrsplit(id,split="__",fixed=TRUE),2),
  #           end=dplyr::nth(data.table::tstrsplit(id,split="__",fixed=TRUE),3),
  #          strand=dplyr::nth(data.table::tstrsplit(id,split="__",fixed=TRUE),4)), by = "id"]
  
  message(" converting to GRanges")
  NF.GR <- GenomicRanges::makeGRangesFromDataFrame(df = NF, keep.extra.columns = TRUE)
  
  message(" cleaning up")
  mcols(NF.GR)[["id"]] <- NULL
  
  message("done!")
  return(NF.GR)
}