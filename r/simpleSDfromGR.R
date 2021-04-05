simpleSDfromGR <- function(
                          ## INPUT
                          GR=NULL,
                          SAMPLE.NAME=NULL,
                          
                          ## OPTIONS
                          USE.READS=TRUE,
                          PLOT.FREQ=TRUE,
                          YLIMS=NULL,
                          RETURN.ALL=FALSE,

                          ## SETTINGS
                          ASPECT.RATIO=1,
                          BAR.FILL="grey80",
                          BAR.LINE="black",
                          XYT=8,
                          FAM="Helvetica",
                          TCOL="black"
                          ){


  ## Pavol Genzor
  ## 11.24.19; Version 1
  ## 02.12.20; Version 2; minor setting adjustments
  ## 03.17.20; Version 3; sequences versus reads, optional output
  
  ## LIBRARIES
  suppressMessages(suppressWarnings(library(data.table)))
  suppressMessages(suppressWarnings(library(ggplot2)))
  suppressMessages(suppressWarnings(library(GenomicRanges)))
  
  ## INPUT CHECKING
  if(is.null(GR)) stop("Please provide a GR  made from .bam file !")
  if(is.null(SAMPLE.NAME)) stop("Please provide SAMPLE.NAME !")
  if(!"MULT" %in% colnames(mcols(GR))) stop("The GRanges must have MULT column!")
  
  ## OPTIONS
  PLOT.ON.Y=ifelse(isTRUE(PLOT.FREQ),"FREQ","COUNT")
  
  ## MAKE TABLE
  GR.DT  <- as.data.table(GR)
  
  ## COUNT SEQ OR READ
  if(isTRUE(USE.READS)){
    message("using reads")
    METRIC.USED = "Read"
    SD.DT <- GR.DT[, lapply(.SD,sum), by = "width", .SDcols = "MULT"] } 
  else {
    message("using sequences")
    METRIC.USED = "Sequence"
    SD.DT <- GR.DT[,.N, by = "width"] }
  
  colnames(SD.DT) <- c("SIZE","COUNT")
  SD.DT[, FREQ := lapply(.SD,function(i){(i/sum(.SD))*100 }),.SDcols = c("COUNT")]
  SD.DT[["SIZE"]] <- as.numeric(SD.DT[["SIZE"]])
  
  ## PLOT
  SD.GG <- ggplot() + theme_pubclean() +
    geom_bar(data = SD.DT, aes_string(x = "SIZE", y = PLOT.ON.Y), 
             stat = "identity", colour = BAR.LINE, fill = BAR.FILL) + 
    scale_x_continuous(breaks = seq(15,50,5)) +
    scale_y_continuous(labels = scientific, 
                       limits = YLIMS, 
                       breaks = seq(0,100,10)) +
    ggtitle(paste(METRIC.USED,SAMPLE.NAME, sep = "; ")) +
    theme(panel.grid = element_blank(), 
          title = element_text(family = FAM, size = XYT, colour = TCOL),
          axis.text = element_text(family = FAM, size = XYT, colour = TCOL),
          axis.title = element_text(family = FAM, size = XYT, colour = TCOL),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          aspect.ratio = ASPECT.RATIO)
  
  ## RETURN
  if(isTRUE(RETURN.ALL)){
    RES.L <- list(SD.DT,SD.GG)
    names(RES.L) <- c(paste0(METRIC.USED,"_data_table"),"ggplot")
    return(RES.L)
  } else { return(SD.GG) }
  
}