
library(packrat)
on()

library(RMassBank)
library(tidyverse)

dir <- "C:/Daten/Michele/20190625 Karin MassBank-QE"
setwd(dir)


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