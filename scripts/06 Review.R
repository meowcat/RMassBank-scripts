library(packrat)
on()
library(here)


source(here("viewer.R"))

w <- loadMsmsWorkspace("results/spectra-mH-eics-score.RData")

review <- viewer(w)

save(review, file="results/review-mH.RData")



w <- loadMsmsWorkspace("results/spectra-pH-eics-score.RData")

review <- viewer(w)

save(review, file="results/review-pH.RData")
