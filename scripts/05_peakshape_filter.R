
source("environment.R")
source("include/functions.R")
library(pROC)
library(RMassBank)
library(gplots)
library(S4Vectors)
library(tidyverse)
library(mzR)
library(glue)

par(mfrow=c(1,1))

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









charge_strs <- adducts

walk(charge_strs, function(charge_str) {
# charge_str <- "pH"

  w <- loadMsmsWorkspace(
    glue("results/spectra-{charge_str}-processed.RData"))



# attach all file names to the compounds
w@spectra <- mapply(function(cpd, f) {
  attr(cpd, "file") <- f
  cpd
}, w@spectra, w@files, SIMPLIFY = FALSE) %>% as("SimpleList")

precursorEps <- 0.1 # maximal m/z deviation of precursor
rtWindow <- 15 # seconds per side
ppmEic <- 5
#ppmGood <- 15

# extract all EICs for MS2 peaks
w@spectra <- lapply(w@spectra, function(cpd) {
  
  selectPolarity <- RMassBank:::getAdductPolarity(charge_str)
  
  f <- attr(cpd, "file")
  message(f)
  d <- openMSfile(f)
  h <- header(d)
  p <- makePeaksCache(d, h)
  
  # 1.6.22: For Exploris data, we need to fill precursor here.
  headerData <- h
  headerData$precursorScanNum <- NA
  headerData[which(headerData$msLevel == 1),"precursorScanNum"] <-
    headerData[which(headerData$msLevel == 1),"acquisitionNum"]
  headerData[,"precursorScanNum"] <- RMassBank:::.locf(headerData[,"precursorScanNum"])
  # Clear the actual MS1 precursor scan number again
  headerData[which(headerData$msLevel == 1),"precursorScanNum"] <- 0
  # Remove precursors which are still NA in precursor scan num.
  # This removes a bug when filling precursor if the first scan(s) are MS2 before a
  # MS1 scan appears. The resulting NA values in precursorScanNum are problematic downstream.
  headerData <- headerData[!is.na(headerData$precursorScanNum),]
  h <- headerData
  
  if(length(cpd@children) == 0)
    return(cpd)
  
  if(!cpd@found)
    return(cpd)
  
  cpd@children <- lapply(cpd@children, function(sp) {
    hSub <- h %>% filter(
      abs(precursorMZ - sp@precursorMz) < precursorEps,
      msLevel == sp@msLevel,
      polarity == selectPolarity,
      abs(retentionTime - sp@rt) < rtWindow,
      collisionEnergy == sp@collisionEnergy
    )
    if(nrow(hSub) == 0)
    {
      attr(sp, "eics") <- list()
      return(sp)
    }
    hSub$msLevel <- 1
      pSub <- p[hSub$seqNum]
    eics <- lapply(property(sp, "mzRaw"), function(mz) {
      eic <- findEIC(d, mz, ppm(mz, ppmEic, p = TRUE), headerCache = hSub, peaksCache = pSub)
      eic$precursorScan <- hSub$precursorScanNum
      eic$mz <- mz
      eic
    })
    
    attr(sp, "eics") <- eics
    return(sp)
  }) %>% as("SimpleList")
  cpd
}) %>% as("SimpleList")

  archiveResults(
    w, 
    glue("results/spectra-{charge_str}-eics.RData")
  )


# correlate EICs with precursor EIC
# note: REMOVE the scan itself and the precursor itself, to zero out the one-point correlations!
w@spectra <- lapply(w@spectra, function(cpd) {
  if(length(cpd@children) == 0)
    return(cpd)
  
  if(!cpd@found)
    return(cpd)

  message(cpd@name)
      
  cpd@children <- lapply(cpd@children, function(sp) {
    eic <- attr(sp, "eic") %>%
      bind_rows(.id = "mzIndex")
    if(nrow(eic) == 0)
    {
      sp <- addProperty(sp, "eicScoreCor", type = "numeric", NA)
      sp <- addProperty(sp, "eicScoreDot", type = "numeric", NA)
      return(sp)
    }
      
    eic <- eic %>%
        dplyr::mutate(mzIndex = as.numeric(as.character(mzIndex))) %>%
        dplyr::select(-mz) %>%
        pivot_wider(names_from = c("mzIndex"), values_from = "intensity", names_prefix = "mz_")
      
      eicPrecursor <- attr(cpd, "eic")
      eic <- eic %>% left_join(eicPrecursor %>% dplyr::rename(precursorScan = scan),
                               by="precursorScan")
      
      
      # filter out "this scan"
      eic <- eic %>% dplyr::filter(precursorScan != sp@precScanNum)
      
      eicPrecursor <- eic[,"intensity"] %>% as.matrix()
      eic <- eic %>% dplyr::select(starts_with("mz_")) %>% as.matrix()
      
    eicScoreCor <- .eicScoreCor(eicPrecursor, eic)
    eicScoreDot <- .eicScoreDot(eicPrecursor, eic)
    
    # # simulate the correlation for a single-point hit
    # singlePointEic <- rep(0, nrow(eicPrecursor))
    # names(singlePointEic) <- rownames(eicPrecursor)
    # singlePointEic[[as.character(sp@precScanNum)]] <- 1
    # singlePointCor <- cor(eicPrecursor, singlePointEic)
    sp <- addProperty(sp, "eicScoreCor", type = "numeric", NA)
    property(sp, "eicScoreCor") <- as.vector(eicScoreCor)
    sp <- addProperty(sp, "eicScoreDot", type = "numeric", NA)
    property(sp, "eicScoreDot") <- as.vector(eicScoreDot)
    # sp@info$singlePointCor <- singlePointCor
    sp
  }) %>% as("SimpleList")
  cpd
}) %>% as("SimpleList")

# evaluation:
# aggregate and keep the best result for every peak
ag <- aggregateSpectra(w)
ag <- ag %>% 
    # group_by(cpdID, scan) %>%
    # # dplyr::mutate(maxint = max(intensity), relint = intensity / max(intensity)) %>%
    # # filter(relint > 0.2) %>%
    # dplyr::mutate(good = abs(dppm) < ppmGood) %>%
  dplyr::group_by(cpdID, scan, mzFound) %>% 
  dplyr::arrange(!good, abs(dppm)) %>% 
  dplyr::slice(1) #%>%
  # dplyr::filter(mz > getOption("RMassBank")$filterSettings$massRangeDivision)
# for every formula, keep the best correlation of all children;
# when no formula present, keep the peak alone
ag$formula_ <- factor(ag$formula) %>% as.numeric()
ag$formula_[is.na(ag$formula_)] <- -seq_along(ag$formula_[is.na(ag$formula_)])
# ag <- ag %>% 
#   group_by(cpdID, formula_) %>% 
#   dplyr::mutate(
#     eicScoreCorMax = max(eicScoreCor, na.rm = TRUE),
#     eicScoreDotMax = max(eicScoreDot, na.rm = TRUE)
#   ) %>%
#   dplyr::mutate(
#     eicScoreCorMin = min(eicScoreCor, na.rm = TRUE),
#     eicScoreDotMin = min(eicScoreDot, na.rm = TRUE)
#   ) 
# 
# ag$eicScoreCorMax[is.infinite(ag$eicScoreCorMax)] <- NA
# ag$eicScoreDotMax[is.infinite(ag$eicScoreDotMax)] <- NA
# ag$eicScoreCorMin[is.infinite(ag$eicScoreCorMin)] <- NA
# ag$eicScoreDotMin[is.infinite(ag$eicScoreDotMin)] <- NA
# # 
# ag$eicScoreCor0 <- ag$eicScoreCor %>% replace_na(0)
# ag$eicScoreDot0 <- ag$eicScoreDot %>% replace_na(0)

ag$eicScoreCor0 <- ag$eicScoreCor %>% replace_na(0)
ag$eicScoreDot0 <- ag$eicScoreDot %>% replace_na(0)

  
#fsc <- f1
fsc <- fBeta(1.5)

setpoint <- roc(good ~ eicScoreCor, data=ag, levels=c("FALSE", "TRUE"), direction="<")
plot(setpoint)
fscore <- fsc(setpoint$sensitivities, setpoint$specificities)
plot(fscore)
maxFscore <- which.max(fscore)
thresholdCor <- setpoint$thresholds[maxFscore]

setpoint <- roc(good ~ eicScoreDot, data=ag, levels=c("FALSE", "TRUE"), direction="<")
plot(setpoint)
fscore <- fsc(setpoint$sensitivities, setpoint$specificities)
  plot(setpoint$thresholds, fscore)
maxFscore <- which.max(fscore)
thresholdDot <- setpoint$thresholds[maxFscore]

setpoint <- roc(good ~ eicScoreCor0, data=ag, levels=c("FALSE", "TRUE"), direction="<")
plot(setpoint)
fscore <- fsc(setpoint$sensitivities, setpoint$specificities)
  plot(setpoint$thresholds, fscore)
maxFscore <- which.max(fscore)
thresholdCor0 <- setpoint$thresholds[maxFscore]


setpoint <- roc(good ~ eicScoreDot0, data=ag, levels=c("FALSE", "TRUE"), direction="<")
plot(setpoint)
fscore <- fsc(setpoint$sensitivities, setpoint$specificities)
  plot(setpoint$thresholds, fscore)
maxFscore <- which.max(fscore)
thresholdDot0 <- setpoint$thresholds[maxFscore]


plot(ag$eicScoreCor, ag$eicScoreDot, col=ag$good+1)
abline(v=thresholdCor, col="red")
abline(h=thresholdDot, col="red")
abline(v=thresholdCor0, col="orange")
abline(h=thresholdDot0, col="orange")





plot(ag$eicScoreCor, ag$eicScoreDot)
hist2d(ag$eicScoreCor, ag$eicScoreDot, FUN = function(x) replace_na(log(length(x)), 0))
hist2d(ag$eicScoreCor, ag$eicScoreDot)

hist(ag$eicScoreCor[ag$good])
hist(ag$eicScoreCor[!ag$good], breaks = 40)
hist(ag$eicScoreCor0[!ag$good], breaks = 40)

hist2d(ag[ag$good,]$eicScoreCor0, ag[ag$good,]$dppm)
abline(v=thresholdCor, col="red")
abline(v=thresholdCor0, col="yellow")
hist2d(ag[ag$good,]$eicScoreCor, log10(ag[ag$good,]$intensity), xlim=c(-1,1))
abline(v=thresholdCor, col="red")
abline(v=thresholdCor0, col="yellow")

hist2d(ag[ag$good,]$eicScoreCor0, log10(ag[ag$good,]$intensity), xlim=c(-1,1))
hist2d(ag[!ag$good,]$eicScoreCor, log10(ag[!ag$good,]$intensity), xlim=c(-1,1))
abline(v=thresholdCor, col="red")
abline(v=thresholdCor0, col="yellow")

hist(ag$eicScoreDot[ag$good])
hist(ag$eicScoreDot[!ag$good], breaks = 40)
hist(ag$eicScoreDot0[!ag$good], breaks = 40)

hist2d(ag[ag$good,]$eicScoreDot0, ag[ag$good,]$dppm)
abline(v=thresholdDot, col="red")
abline(v=thresholdDot0, col="yellow")
hist2d(ag[ag$good,]$eicScoreDot, log10(ag[ag$good,]$intensity), xlim=c(0,1))
abline(v=thresholdDot, col="red")
abline(v=thresholdDot0, col="yellow")

hist2d(ag[ag$good,]$eicScoreDot0, log10(ag[ag$good,]$intensity), xlim=c(0,1))
abline(v=thresholdDot, col="red")
abline(v=thresholdDot0, col="yellow")
hist2d(ag[!ag$good,]$eicScoreDot, log10(ag[!ag$good,]$intensity), xlim=c(0,1))
abline(v=thresholdDot, col="red")
abline(v=thresholdDot0, col="yellow")



plot(ag[!ag$good,]$eicScoreCor, ag[!ag$good,]$eicScoreDot)
plot(ag[ag$good,]$eicScoreCor, ag[ag$good,]$eicScoreDot)

hist2d(ag[!ag$good,]$eicScoreCor, ag[!ag$good,]$eicScoreDot)
abline(v=thresholdCor, col="red")
abline(h=thresholdDot, col="red")
abline(v=thresholdCor0, col="orange")
abline(h=thresholdDot0, col="orange")

hist2d(ag[ag$good,]$eicScoreCor, ag[ag$good,]$eicScoreDot, 
       FUN=length)
abline(v=thresholdCor, col="red")
abline(h=thresholdDot, col="red")
abline(v=thresholdCor0, col="orange")
abline(h=thresholdDot0, col="orange")


plot(ag[!ag$good,]$eicScoreCor, ag[!ag$good,]$eicScoreDot)
abline(v=thresholdCor, col="red")
abline(h=thresholdDot, col="red")
abline(v=thresholdCor0, col="orange")
abline(h=thresholdDot0, col="orange")

plot(ag[ag$good,]$eicScoreCor, ag[ag$good,]$eicScoreDot)
abline(v=thresholdCor, col="red")
abline(h=thresholdDot, col="red")
abline(v=thresholdCor0, col="orange")
abline(h=thresholdDot0, col="orange")

attr(w, "eicScoreFilter") <- list("eicScoreCor" = thresholdCor,
                                 "eicScoreDot" = thresholdDot)

  archiveResults(
    w, 
    glue("results/spectra-{charge_str}-eics-score.RData")
  )
  
   }) # walk
