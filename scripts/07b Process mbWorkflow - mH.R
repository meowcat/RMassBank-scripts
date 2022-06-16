library(here)

library(RMassBank)
library(tidyverse)

loadList("results/compoundsRmb.csv")
loadRmbSettings("input/RmbSettings.ini")


#w <- loadMsmsWorkspace("results/spectra-pH-processed.RData")

load("results/review-mH.RData")

# adjust the data in the spectra as long as there is no other way to select it
w_ <- review$w
for(i in seq_along(w_@spectra))
  attr(w_@spectra[[i]], "specOK") <- review$specOk[[i]]
w_@spectra <- w_@spectra[review$cpdOk]

metric <- "eicScoreCor"
eicScoreLimit <- attr(w_, "eicScoreFilter")[[metric]]

w_@spectra <- as(lapply(w_@spectra, function(cpd) {
  specOK <- attr(w_, "specOK")
  for(i in seq_along(specOK))
    cpd@children[[i]]@ok <- specOK[[i]]
  cpd@children <- as(lapply(cpd@children, function(sp) {
    d <- getData(sp)
    d <- d %>% group_by(mz) %>% arrange(!good, abs(dppm)) %>% slice(1)
    d <- d %>% mutate(
      good_ = good,
      good = eicScoreCor > eicScoreLimit
    ) %>% filter(good)
    sp <- setData(sp, as.data.frame(d))
    sp
  }), "SimpleList")
  cpd
}), "SimpleList")

mb <- newMbWorkspace(w_)
mb <- resetInfolists(mb)
mb <- loadInfolist(mb, "infolist_checked.csv")

mb <- mbWorkflow(mb, steps=c(1:3), infolist_path = "infolist_mH.csv")
save(mb, file="results/records_mH.RData")

