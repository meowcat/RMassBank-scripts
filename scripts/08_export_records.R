library(here)

library(RMassBank)
library(tidyverse)
library(glue)

RMassBank.env$export.molfiles <- FALSE


loadList("results/compoundsRmb.csv")
loadRmbSettings("input/RmbSettings.ini")

fs::dir_create("results/records")
fs::dir_delete("results/records")
fs::dir_create("results/records")

charge_strs <- c("mH", "pH")
walk(charge_strs, function(charge_str) {
  
  load(glue("results/records_{charge_str}.RData"))
  mb <- resetInfolists(mb)
  mb <- loadInfolists(mb, "infolists")
  
  wd <- getwd()
  setwd("results/records")
  mb <- mbWorkflow(mb)
  setwd(wd)
  
})
