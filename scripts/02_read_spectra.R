
library(RMassBank)
library(mzR)
library(plyr)
library(dplyr)

source("environment.R")

if(!exists("column_of_id"))
  column_of_id <- 3

inFiles <- list.files(inputDir, ".mzML", full.names = T)

compounds <- read.csv(compoundsCsv, sep=",")

checkFiles <- data.frame(files = inFiles, stringsAsFactors = FALSE)
checkFiles <- checkFiles %>% mutate(basefiles = basename(files)) %>%
  mutate(id = strsplit(basefiles, '[_.]') %>% lapply(`[`, column_of_id) %>% as.integer)

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

charge_strs <- adducts

walk(charge_strs, function(charge_str) {
  
  
  
  w <- newMsmsWorkspace()
  w <- msmsRead(w, filetable="results/compoundsMsR.csv", 
                readMethod="mzR", mode=charge_str )
  eics <- alply(compoundsMsRead, 1, function(cpd)
  {
    message(cpd[["Files"]])
    d <- openMSfile(cpd[["Files"]])
    h <- header(d)
    h <- h[h$polarity == RMassBank:::getAdductPolarity(charge_str),]
    mz <- findMz(cpd[["ID"]], mode = charge_str, ppm = 5)
    findEIC(d, mz, headerCache = h)
  })
  
  
  for(i in seq_along(w@spectra))
  {
    attr(w@spectra[[i]], "eic") <- eics[[i]]
  }
  archiveResults(w, glue::glue(
    "results/spectra-{charge_str}-msmsRead.RData"
    ))
  
})


