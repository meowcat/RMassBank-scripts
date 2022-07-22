library(here)

library(RMassBank)
library(tidyverse)
library(glue)
library(fs)
library(assertthat)

source(here::here("functions.R"))

loadList("results/compoundsRmb.csv")
loadRmbSettings("input/RmbSettings.ini")



infolists <- c()

charge_strs <- c("pH", "mH")

# Check that the necessary files exist
walk(charge_strs, function(charge_str) {
  assert_that(fs::file_exists(glue("results/review_{charge_str}_score_cutoff.csv")),
              msg = glue("You have not yet reviewed {charge_str} data")
  )
  assert_that(fs::file_exists(glue("results/review_{charge_str}_cpd_ok.csv")),
              msg = glue("You have not yet reviewed {charge_str} data")
  )
  assert_that(fs::file_exists(glue("results/review_{charge_str}_spec_ok.csv")),
              msg = glue("You have not yet reviewed {charge_str} data")
  )
})


# run mbworkflow (collect online data) for all specified modes
walk(charge_strs, function(charge_str) {
  
  
  cpdOk <- read_csv(glue("results/review_{charge_str}_cpd_ok.csv"))
  specOk <- read_csv(glue("results/review_{charge_str}_spec_ok.csv")) # %>%
    # group_by(cpd) %>%
    # group_split() %>%
    # map(~.x$ok)
  
  cutoff_file <- glue("results/review_{charge_str}_score_cutoff.csv")
  if(fs::file_exists(cutoff_file))
    score_cutoff <- read_file(cutoff_file) %>% as.numeric()
  else
    score_cutoff <- NA
  
  
  w_ <- loadMsmsWorkspace(glue("results/spectra-{charge_str}-autoreview.RData"))
  
  # Inject review data into compounds,
  # then select and deselect spectra based on the choices from the review files
  for(i in seq_along(w_@spectra)) {
    attr(w_@spectra[[i]], "specOK") <- specOk %>% filter(cpd == i)
    # assert_that(all(attr(w_@spectra[[i]], "specOK")$name == ))
    attr(w_@spectra[[i]], "threshold") <- cpdOk$threshold[[i]]
  }
  w_@spectra <- w_@spectra[cpdOk$ok]
  w_@files <- w_@files[cpdOk$ok]
  
  # Have a fallback in case eicScoreCor was not set for whatever reason
  metric <- "eicScoreCor"
  if(!is.na(score_cutoff)) {
    eicScoreLimit <- score_cutoff
  } else {
    warning("score cutoff could not be read from review data, using calculated default")
    eicScoreLimit <- attr(w_, "eicScoreFilter")[[metric]]
  }
    
  
  w_@spectra <- as(lapply(w_@spectra, function(cpd) {
    specOK <- attr(cpd, "specOK")
    assert_that(nrow(specOK) == length(cpd@children))
  
    for(i in seq_along(cpd@children)) {
      cpd@children[[i]]@ok <- specOK$ok[[i]]
      attr(cpd@children[[i]], "threshold") <- specOK$threshold[[i]]
    }
    cpd@children <- as(lapply(cpd@children, function(sp) {
      
      # find the valid threshold for this spectrum:
      # can be the global, compound-specific or spectrum-specific threshold
      thresholds <- c(
        eicScoreLimit,
        attr(cpd, "threshold"),
        attr(sp, "threshold")
      )
      threshold <- tail(thresholds[!is.na(thresholds)], 1)
      if(!sp@ok)
        return(sp)
      
      d <- getData(sp)
      d <- d %>% group_by(mz) %>% arrange(!good, desc(formulaMultiplicity), abs(dppm)) %>% slice(1)
      d <- d %>% mutate(
        good_ = good,
        good = (eicScoreCor > threshold) & !satellite & !low
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