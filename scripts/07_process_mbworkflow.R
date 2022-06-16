library(here)

library(RMassBank)
library(tidyverse)
library(glue)

source(here("functions.R"))

loadList("results/compoundsRmb.csv")
loadRmbSettings("input/RmbSettings.ini")



infolists <- c()

charge_strs <- c("pH", "mH")

walk(charge_strs, function(charge_str) {
  
  
  cpdOk <- read_csv(glue("input/review_{charge_str}_cpd_ok.csv")) %>% pull(ok)
  specOk <- read_csv(glue("input/review_{charge_str}_spec_ok.csv")) %>%
    group_by(cpd) %>%
    group_split() %>%
    map(~.x$ok)
  
  cutoff_file <- glue("input/review_{charge_str}_score_cutoff.csv")
  if(fs::file_exists(cutoff_file))
    score_cutoff <- read_file(cutoff_file) %>% as.numeric()
  else
    score_cutoff <- NA
  
  
  w_ <- loadMsmsWorkspace(glue("results/spectra-{charge_str}-autoreview.RData"))
  
  # Select and deselect spectra based on the choices from the review files
  for(i in seq_along(w_@spectra))
    attr(w_@spectra[[i]], "specOK") <- specOk[[i]]
  w_@spectra <- w_@spectra[cpdOk]
  
  metric <- "eicScoreCor"
  if(!is.na(score_cutoff))
    eicScoreLimit <- score_cutoff
  else
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
  for(infolist in infolists)
    mb <- loadInfolist(mb, infolist)
  newInfolist <- glue("results/infolist_{charge_str}.csv")
  mb <- mbWorkflow(mb, steps=c(1:3), 
                   infolist_path = newInfolist)
  
  infolists <<- c(infolists, newInfolist)
  save(mb, file=glue("results/records_{charge_str}.RData"))
})