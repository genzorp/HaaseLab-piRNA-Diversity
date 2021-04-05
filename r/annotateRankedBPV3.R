annotateRankedBPV3 <- function(
                                ## INPUT
                                ANN.TAB.L = NULL,

                                ## SETTINGS
                                ASPECT.RATIO = 2,
                                Y.LAB = "Fraction (%)",
                                FAM = "Helvetica", TCOL = "black", XYT = 12,
                                RETURN.ALL = FALSE){
  
  ## Pavol Genzor
  ## To plot sense and antisense annotations
  ## 05.01.20; Version 1
  
  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("dplyr"); library("ggplot2")})
  
  ## INPUT
  if(is.null(ANN.TAB.L)) stop("Please provide ANN.TAB.L !")
  if(!is.list(ANN.TAB.L)) stop("Please provide list of tables !")
  
  ## Combine samples into directional tables
  DT.UP <- rbindlist(lapply(ANN.TAB.L, function(S){S[["DTM.UP"]]}))
  DT.UP[["SAMPLE"]] <- factor(DT.UP[["SAMPLE"]], levels = unique(DT.UP[["SAMPLE"]]))
  DT.DOWN <- rbindlist(lapply(ANN.TAB.L, function(S){S[["DTM.DOWN"]]}))
  DT.DOWN[["SAMPLE"]] <- factor(DT.DOWN[["SAMPLE"]], levels = unique(DT.DOWN[["SAMPLE"]]))
  
  ## All data table
  DT <- rbindlist(list(DT.UP,DT.DOWN))
  
  ## Plot
  GGBP <- ggplot() + theme_bw() +
    ## DATA
    geom_bar(data = DT.UP, aes(x = SAMPLE, y = UD.VALUE), 
             stat = "identity", fill = DT.UP[["COLS"]], colour = NA) + 
    geom_bar(data = DT.DOWN, aes(x = SAMPLE, y = UD.VALUE), 
             stat = "identity", fill = DT.DOWN[["COLS"]], colour = NA) + 
    
    ## LABELS
    geom_text(data = DT.UP, aes(x = SAMPLE, y = UD.VALUE, label = NAME), 
              position = position_stack(vjust = 0.5), 
              colour=TCOL, family=FAM, size = 2) +
    geom_text(data = DT.DOWN[NAME %in% "OTHER"], aes(x = SAMPLE, y = UD.VALUE, label = NAME), 
              position = position_stack(vjust = 0.5), 
              colour="white", family=FAM, size = 2) +
    
    ## LABS
    xlab("") + ylab(paste0("Fraction (%)")) + 
    ggtitle(paste0("Annotation Barplot\n",
                   "used: ",unique(DT[["dataType"]]),"\n",
                   "filters: ",unique(DT[["FILTERS"]]))) +
    
    ## LINE
    geom_hline(yintercept = 0, size = 0.5, colour = "black" ) +
    
    ## SCALES
    scale_y_continuous(limits = c(-100,100), breaks = seq(-100,100,10)) +
    
    ## THEME
    theme(aspect.ratio = ASPECT.RATIO, 
          panel.border = element_blank(), panel.grid = element_blank(),
          axis.text.x = element_text(family = FAM, color = TCOL, size = XYT, 
                                     angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(family = FAM, color = TCOL, size = XYT),
          axis.title = element_text(family = FAM, color = TCOL, size = XYT))
  
  
  ## RESULTS
  RES.L <- list("plotData" = DT, "ggplot" = GGBP)
  
  ## RETURN
  if(isTRUE(RETURN.ALL)) {
    return(RES.L) } else { return(RES.L[["ggplot"]]) }
  
}