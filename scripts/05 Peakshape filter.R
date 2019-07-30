

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

# Peak EIC correlation while filtering out zero rows
.eicScoreCor <- function(prec, eic) {
  eicSum <- colSums(abs(eic))
  .eicScore <- cor(prec, eic[,eicSum > 0,drop=FALSE])
  eicScore <- rep(NA_real_, ncol(eic))
  eicScore[eicSum > 0] <- .eicScore
  eicScore
} 

# Peak EIC dot product score
.eicScoreDot <- function(prec, eic) {
  normX <- sqrt(sum(prec^2))
  normY <- sqrt(colSums(eic^2))
  dot <- t(prec) %*% eic
  score <- dot / (normX * normY)
  score[normY == 0] <- NA
  score
}


# correlate EICs with precursor EIC
# note: REMOVE the scan itself and the precursor itself, to zero out the one-point correlations!
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
    # assert that the rows match
    if(!all(rownames(eic) == rownames(eicPrecursor)))
      stop("EIC correlations: the scan index of precursor and spectrum EIC doesn't match!")
    # remove actual scan and precursor
    # don't run columns that have zero intensity outside of the self-point,
    # and replace with NA
    thisScan <- which(rownames(eic) == as.character(sp@precScanNum))
    eic <- eic[-thisScan,,drop=FALSE]
    eicPrecursor <- eicPrecursor[-thisScan,,drop=FALSE]
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
# for every formula, keep the best correlation of all children
ag <- aggregateSpectra(w)
ag <- ag %>% group_by(cpdID, scan, mzFound) %>% arrange(!good, abs(dppm)) %>% slice(1)
ag <- ag %>% ungroup() %>% group_by(cpdID, formula) %>% mutate(
  eicScoreCor = max(eicScoreCor, na.rm = TRUE), eicScoreDot = max(eicScoreDot, na.rm = TRUE)
) 

ag$eicScoreCor0 <- ag$eicScoreCor %>% replace_na(0)
ag$eicScoreDot0 <- ag$eicScoreDot %>% replace_na(0)


setpoint <- roc(good ~ eicScoreCor, data=ag, levels=c("FALSE", "TRUE"), direction="<")
plot(setpoint)
fscore <- 2 * setpoint$sensitivities * setpoint$specificities / (setpoint$sensitivities + setpoint$specificities)
plot(fscore)
maxFscore <- which.max(fscore)
thresholdCor <- setpoint$thresholds[maxFscore]

setpoint <- roc(good ~ eicScoreDot, data=ag, levels=c("FALSE", "TRUE"), direction="<")
plot(setpoint)
fscore <- 2 * setpoint$sensitivities * setpoint$specificities / (setpoint$sensitivities + setpoint$specificities)
plot(fscore)
maxFscore <- which.max(fscore)
thresholdDot <- setpoint$thresholds[maxFscore]

ag$eicScoreCor0 <- ag$eicScoreCor %>% replace_na(0)
ag$eicScoreDot0 <- ag$eicScoreDot %>% replace_na(0)


setpoint <- roc(good ~ eicScoreCor0, data=ag, levels=c("FALSE", "TRUE"), direction="<")
plot(setpoint)
fscore <- 2 * setpoint$sensitivities * setpoint$specificities / (setpoint$sensitivities + setpoint$specificities)
plot(fscore)
maxFscore <- which.max(fscore)
thresholdCor0 <- setpoint$thresholds[maxFscore]


setpoint <- roc(good ~ eicScoreDot0, data=ag, levels=c("FALSE", "TRUE"), direction="<")
plot(setpoint)
fscore <- 2 * setpoint$sensitivities * setpoint$specificities / (setpoint$sensitivities + setpoint$specificities)
plot(fscore)
maxFscore <- which.max(fscore)
thresholdDot0 <- setpoint$thresholds[maxFscore]


plot(ag$eicScoreCor, ag$eicScoreDot)
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

.lenWrap <- function(fun) function(x) {
  res <- fun(x)
  if(length(res) == 0) res <- 0
  if(is.na(res)) res <- 0
  res
}

hist2d(ag$eicScoreCor, ag$eicScoreDot, FUN = .lenWrap(function(x) log(length(x))))




points(ag[!ag$good,]$eicScoreCor, ag[!ag$good,]$eicScoreDot, col="red")
hist2d(ag$eicScoreCor, ag$eicScoreDot, FUN = function(x) replace_na(log(length(x)), 0))
hist2d(ag$eicScoreCor, ag$eicScoreDot)

# not useful, since there are very few of these peaks
# hist2d(ag[!ag$good,]$eicScore0, ag[!ag$good,]$dppm)
# not useful, since almost all of them are 0 then_
# hist2d(ag[!ag$good,]$eicScore, replace_na(ag[!ag$good,]$dppm, 0))

# w <- loadMsmsWorkspace("results/spectra-mH-msmsRead.RData")
# w <- msmsWorkflow(w, "mH", c(2:4))
# archiveResults(w, "results/spectra-mH-processed-step4.RData")
# w <- msmsWorkflow(w, "mH", c(5:8))
# archiveResults(w, "results/spectra-mH-processed.RData")
