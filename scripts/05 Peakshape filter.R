

library(RMassBank)
library(reshape2)

dir <- "C:/Daten/Michele/20190625 Karin MassBank-QE/"
setwd(dir)


loadList("results/compoundsRmb.csv")
loadRmbSettings("input/RmbSettings.ini")

# 
# recalibrate.lumos <- function (rcdata) 
# {
#   span <- 0.25
#   mingroups <- nrow(rcdata[!is.na(rcdata$mzFound), ])
#   if (mingroups < 4) {
#     warning("recalibrate.loess: Not enough data points, omitting recalibration")
#     return(recalibrate.identity(rcdata))
#   }
#   else if (mingroups * span < 4) {
#     span <- 4/mingroups
#     warning("recalibrate.loess: Span too small, resetting to ", 
#             round(span, 2))
#   }
#   return(loess(recalfield ~ mzFound, data = rcdata, family = c("symmetric"), 
#                degree = 1, span = 0.4, surface = "direct"))
# }




w <- loadMsmsWorkspace("results/spectra-pH-processed.RData")


library(S4Vectors)
library(tidyverse)
library(mzR)

# attach all file names to the compounds
w@spectra <- mapply(function(cpd, f) {
  attr(cpd, "file") <- f
  cpd
}, w@spectra, w@files, SIMPLIFY = FALSE) %>% as("SimpleList")

precursorEps <- 0.1 # maximal m/z deviation of precursor
rtWindow <- 15 # seconds per side
ppmEic <- 5
selectPolarity <- 1

# extract all EICs for MS2 peaks
w@spectra <- lapply(w@spectra, function(cpd) {
  f <- attr(cpd, "file")
  message(f)
  d <- openMSfile(f)
  h <- header(d)
  p <- makePeaksCache(d, h)
  if(length(cpd@children) == 0)
    return(cpd)
  cpd@children <- lapply(cpd@children, function(sp) {
    hSub <- h %>% filter(
      abs(precursorMZ - sp@precursorMz) < precursorEps,
      msLevel == sp@msLevel,
      polarity == selectPolarity,
      abs(retentionTime - sp@rt) < rtWindow,
      collisionEnergy == sp@collisionEnergy
    )
    hSub$msLevel <- 1
    pSub <- p[hSub$acquisitionNum]
    eics <- lapply(property(sp, "mzRaw"), function(mz) {
      eic <- findEIC(d, mz, ppm(mz, ppmEic, p = TRUE), headerCache = hSub, peaksCache = pSub)
      eic$precursorScan <- hSub$precursorScanNum
      eic
    })
    
    attr(sp, "eics") <- eics
    return(sp)
  }) %>% as("SimpleList")
  cpd
}) %>% as("SimpleList")

archiveResults(w, "results/spectra-pH-eics.RData")

# correlate EICs with precursor EIC
w@spectra <- lapply(w@spectra, function(cpd) {
  if(length(cpd@children) == 0)
    return(cpd)
  cpd@children <- lapply(cpd@children, function(sp) {
    eic <- attr(sp, "eic") %>%
      bind_rows(.id = "mzIndex") %>%
      acast(precursorScan ~ as.numeric(as.character(mzIndex)), value.var = "intensity",
            drop = FALSE)
    eicPrecursor <- attr(cpd, "eic") %>% 
      filter(scan %in% rownames(eic)) %>%
      acast(scan ~ 1, value.var = "intensity")
    eicCor <- cor(eicPrecursor, eic)
    # simulate the correlation for a single-point hit
    singlePointEic <- rep(0, nrow(eicPrecursor))
    names(singlePointEic) <- rownames(eicPrecursor)
    singlePointEic[[as.character(sp@precScanNum)]] <- 1
    singlePointCor <- cor(eicPrecursor, singlePointEic)
    sp <- addProperty(sp, "cor", type = "numeric", NA)
    property(sp, "cor") <- as.vector(eicCor)
    sp@info$singlePointCor <- singlePointCor
    sp
  }) %>% as("SimpleList")
  cpd
}) %>% as("SimpleList")

# evaluation:
# Mode of unmatched peaks



# w <- loadMsmsWorkspace("results/spectra-mH-msmsRead.RData")
# w <- msmsWorkflow(w, "mH", c(2:4))
# archiveResults(w, "results/spectra-mH-processed-step4.RData")
# w <- msmsWorkflow(w, "mH", c(5:8))
# archiveResults(w, "results/spectra-mH-processed.RData")
