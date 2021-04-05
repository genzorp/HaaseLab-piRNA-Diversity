simpleSRViolin_V5 <- function(## INPUT
                              GRL = NULL,
                              
                              ## OPTIONS
                              Y.PPM = FALSE,
                              
                              ## SETTINGS
                              SAMPLE.ORDER = names(GRL),
                              SIZE.RANGE = c(18,50), 
                              NH.TAG = NULL,
                              
                              ## OUTPUT
                              ADD.TABLE = TRUE,
                              SIG.FIGURES = 2,
                              SOURCE.DIR = NULL,
                              RETURN.ALL = FALSE,
                              
                              ## PLOT SETTINGS
                              Y.LIMS = NULL,
                              VIOLIN.ADJUST = 1, 
                              VIOLIN.SCALE = "width",
                              ASPECT.RATIO = 1,
                              LEGEND.POSITION = "right",
                              TCOL = "black",FAM = "Helvetica",XYT = 12,
                              TRANS = "log10",PLOT.COLORS = c("orange1","grey70")){
  
  ## Pavol Genzor
  ## Based on scripts from Daniel Stoyko
  ## 03.18.20; Version 1
  ## 04.01.20; Version 2; Adding table
  ## 04.03.20; Version 3; Adding y.limits
  ## 05.04.20; Version 4; added consistent filter and library loading
  ## 07.06.20; Version 5; adding PPM conversion (some of Daniel's edits)
  
  ## LIBRARIES
  suppressPackageStartupMessages({
    library("data.table"); library("dplyr"); library("GenomicRanges")
    library("ggplot2"); library("tidyverse"); library("reshape")})
  
  ## INPUT CHECKING
  if(is.null(GRL)) stop("Please provide a named GRL object !")
  if(isFALSE(is.list(GRL))) stop("Input is not a GRL - Please proide a list of GRanges !")
  if(is.null(names(GRL))) stop("Please make sure that GRL is named !")
  if(is.null(SOURCE.DIR)) stop("Please provide a SOURCE.DIR with functions !")
  if(!"MULT" %in% colnames(mcols(GRL[[1]]))) stop("GRL needs to have MULT column !")
  if(!"NH" %in% colnames(mcols(GRL[[1]]))) stop("GRL needs to have NH column !")
  
  ## FUNCTION
  source(paste0(SOURCE.DIR,"simpleGRFilter.R"))
  
  ## PROCEED
  message("Processing...")

  ##
  ## DATA LOOP
  ##
  
  GRL.DT.L <- lapply(names(GRL), function(i){
    
    ## SUBSET
    GR <- GRL[[i]]
    
    ## FILTER
    GR <- simpleGRFilter(GR = GR, 
                         RANGE.NAME = i,
                         NH.TAG = NH.TAG, 
                         SIZE.RANGE = SIZE.RANGE)
    
    ## PREP THE TABLE
    GR.DT <- as.data.table(GR)
    NDT <- GR.DT[,.N, by = "MULT"]
    setorderv(x = NDT, cols = "MULT")
    colnames(NDT) <- c("MULT","SEQ")
    
    ## CALCULATE 
    NDT[, READ := MULT * SEQ]
    NDT[, SAMPLE := i]

    ## REARRANGE
    mNDT <- setDT(melt(data = NDT,
                       id.vars = c("MULT","SAMPLE"), 
                       variable.name = "TYPE", 
                       value.name = "COUNT"))
    ## RETURN
    return(mNDT) })
  
  ## COMNBINE AND ORGANIZE FOR PLOTTING
  message(" calculating & plotting")
  GR.DT <- rbindlist(GRL.DT.L)
  GR.DT[["SAMPLE"]] <- factor(GR.DT[["SAMPLE"]], levels = SAMPLE.ORDER)
  GR.DT[["TYPE"]] <- factor(GR.DT[["TYPE"]], levels = c("SEQ","READ"))
  GR.DT[["GROUP"]] <- paste(GR.DT[["SAMPLE"]],GR.DT[["TYPE"]], sep = "_")
  for(i in unique(GR.DT[["GROUP"]])) set(x = setDT(GR.DT), i = which(GR.DT[["GROUP"]] %in% i), 
                                         j = "GROUP_SUM", value = sum(GR.DT[GROUP %in% i][["COUNT"]]))
  
  ##
  ## TABLE FOR PLOT
  ##
  
  TAB.DT <- GR.DT[, sum(.SD), by = c("SAMPLE","TYPE"), .SDcols = "COUNT"]
  TAB.DT[["POOL"]] <- nth(tstrsplit(TAB.DT[["SAMPLE"]], split = "-"),-1)
  for(t in unique(TAB.DT[["TYPE"]])) set(x = setDT(TAB.DT), 
                                         i = which(TAB.DT[["TYPE"]] %in% t), 
                                         j = "FRACTION", 
                                         value = round(TAB.DT[TYPE %in% t][["V1"]]/
                                                         TAB.DT[TYPE %in% t & POOL %in% "total"][["V1"]],
                                                       digits = SIG.FIGURES))
  TAB.W <- dcast.data.table(data = TAB.DT, formula = TYPE ~ SAMPLE, value.var = "FRACTION")

  ## CONVERT TO PPM
  if(isTRUE(Y.PPM)){
    message("... MULT to PPM conversion (adjust Y.LIMS) ")
    GR.DT[["MULT"]] <- (GR.DT[["MULT"]]/TAB.DT[TYPE %in% "READ" & POOL %in% "total"][["V1"]])*1000000
    Y.LAB <- ylab("MULT.PPM")
    } else { Y.LAB <- ylab("MULT") }
  
  ##
  ## VIOLIN PLOT
  ##
  
  GG.VIOL <-  ggplot() + theme_pubclean() +
    ## DATA
    geom_violin(data = GR.DT, aes(x = SAMPLE, y = MULT, weight = COUNT/GROUP_SUM, fill = TYPE),
                color="black", lwd = 0.25, adjust= VIOLIN.ADJUST, scale = VIOLIN.SCALE, 
                position = position_dodge(width = 1)) +
    
    ## TITLE
    ggtitle(paste0("Violins settings\n",
                   paste0("size range: ", min(SIZE.RANGE)," - ", max(SIZE.RANGE)),"\n",
                   paste0("used NH tags: ", ifelse(is.null(NH.TAG),"all",NH.TAG) ))) + Y.LAB +
    
    ## SCALES
    scale_fill_manual(values = PLOT.COLORS) +
    scale_y_continuous(trans = TRANS, limits = Y.LIMS) +
    annotation_logticks(sides = "l") +
    theme(aspect.ratio = ASPECT.RATIO, 
          axis.ticks.y = element_blank(),
          legend.position = LEGEND.POSITION,
          panel.grid = element_blank(), 
          axis.text = element_text(family = FAM, color =TCOL, size = XYT))
  
  ## CONVERT INTO GROP
  if(isTRUE(ADD.TABLE)){
    TBL <- tableGrob(d = TAB.W, theme = ttheme_minimal(), rows = NULL)
    GG.VIOL <- grid.arrange(GG.VIOL,TBL, nrow = 2, heights = c(5,1))}
  
  ## RETURN
  message("Done.")
  if(isTRUE(RETURN.ALL)){
    RES.L <- list("plotData" = GR.DT,
                  "ggplot" = GG.VIOL)
    return(RES.L)
  } else { return(GG.VIOL) }
}