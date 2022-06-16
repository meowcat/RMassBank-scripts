library(here)

source("environment.R")

library(RMassBank)

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


o <- getOption("RMassBank")
o$electronicNoise <- c()
options(RMassBank = o)

w <- loadMsmsWorkspace("results/spectra-pH-msmsRead.RData")
w <- msmsWorkflow(w, "pH", c(2:4))
archiveResults(w, "results/spectra-pH-processed-step4.RData")
w <- msmsWorkflow(w, "pH", c(5:8))
archiveResults(w, "results/spectra-pH-processed.RData")

w <- loadMsmsWorkspace("results/spectra-mH-msmsRead.RData")
w <- msmsWorkflow(w, "mH", c(2:4))
archiveResults(w, "results/spectra-mH-processed-step4.RData")
w <- msmsWorkflow(w, "mH", c(5:8))
archiveResults(w, "results/spectra-mH-processed.RData")
