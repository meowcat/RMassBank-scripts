
library(here)


library(RMassBank)
library(mzR)
library(plyr)
library(dplyr)
# dir <- "C:/Daten/Michele/20190625 Karin MassBank-QE/"
# setwd(dir)

source("environment.R")

inFiles <- list.files(inputDir, ".mzML", full.names = T)

compounds <- read.csv(compoundsCsv, sep=",")

checkFiles <- data.frame(files = inFiles, stringsAsFactors = FALSE)
checkFiles <- checkFiles %>% mutate(basefiles = basename(files)) %>%
  mutate(id = strsplit(basefiles, '[_.]') %>% lapply(`[`, 3) %>% as.integer)

compounds$found <- compounds$ID %in% checkFiles$id



compoundsRmb <- compounds[compounds$found,
                          c("ID", "Compound", "SMILES", "RT")]
#compoundsRmb$RT <- ""
compoundsRmb$CAS <- ""
colnames(compoundsRmb) <- c("ID", "Name", "SMILES", "RT", "CAS")
write.csv(compoundsRmb, file="results/compoundsRmb.csv")

loadList("results/compoundsRmb.csv")
loadRmbSettings("input/RmbSettings.ini")

compounds <- merge(compounds, checkFiles, by.x='ID', by.y='id', all.x = TRUE)

compoundsMsRead <- compounds[compounds$found,c("ID", "files")]
#compoundsMsRead$file <- paste(inputDir, compoundsMsRead$file, sep="/")
colnames(compoundsMsRead) <- c("ID", "Files")
write.csv(compoundsMsRead, file="results/compoundsMsR.csv")

# Process positive mode
w <- newMsmsWorkspace()
#wPeak <- msmsRead.peak(w, filetable="results/compoundsMsR.csv", readMethod="mzR", mode="pH" )
w <- msmsRead(w, filetable="results/compoundsMsR.csv", readMethod="mzR", mode="pH" )

eics <- alply(compoundsMsRead, 1, function(cpd)
  {
  message(cpd[["Files"]])
  d <- openMSfile(cpd[["Files"]])
  h <- header(d)
  h <- h[h$polarity == 1,]
  mz <- findMz(cpd[["ID"]], mode = "pH", ppm = 5)
  findEIC(d, mz, headerCache = h)
})

for(i in seq_along(w@spectra))
{
  attr(w@spectra[[i]], "eic") <- eics[[i]]
}
archiveResults(w, "results/spectra-pH-msmsRead.RData")

# Process negative mode

w <- newMsmsWorkspace()
#wPeak <- msmsRead.peak(w, filetable="results/compoundsMsR.csv", readMethod="mzR", mode="mH" )
w <- msmsRead(w, filetable="results/compoundsMsR.csv", readMethod="mzR", mode="mH" )
eics <- alply(compoundsMsRead, 1, function(cpd)
{
  message(cpd[["Files"]])
  d <- openMSfile(cpd[["Files"]])
  h <- header(d)
  h <- h[h$polarity == 0,]
  mz <- findMz(cpd[["ID"]], mode = "mH", ppm = 10)
  findEIC(d, mz, headerCache = h)
})

for(i in seq_along(w@spectra))
{
  attr(w@spectra[[i]], "eic") <- eics[[i]]
}
archiveResults(w, "results/spectra-mH-msmsRead.RData")




