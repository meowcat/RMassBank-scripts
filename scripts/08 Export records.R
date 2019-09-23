library(packrat)
on()
library(here)

library(RMassBank)
library(tidyverse)

loadList("results/compoundsRmb.csv")
loadRmbSettings("input/RmbSettings.ini")


load("results/records_pH.RData")
mb <- resetInfolists(mb)
mb <- loadInfolist(mb, "infolist_checked.csv")
mb <- mbWorkflow(mb, filter = FALSE)


load("results/records_mH.RData")
mb <- resetInfolists(mb)
mb <- loadInfolist(mb, "infolist_checked.csv")
mb <- mbWorkflow(mb, filter = FALSE)