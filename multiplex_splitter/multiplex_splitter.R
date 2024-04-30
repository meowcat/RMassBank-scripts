source("environment.R")


library(rawrr)
library(tidyverse)
library(fs)
library(glue)

rawfiles <- fs::dir_ls(inputMultiplex, glob = "*.mzML")
mapping <- read_csv(fs::path(inputMultiplex, "mapping.csv"), col_names = c("raw_file", "compound_map", "batch"))

fs::dir_create("output/compounds_map")
fs::dir_create(inputDir, recurse = TRUE)

cosomi_data <- read_csv(fs::path(inputDir, "compounds_cosomi.csv"))

# Match compounds to raw files
tables <- pmap_dfr(mapping, function(raw_file, compound_map, batch) {

  map_data <- read_csv(fs::path(inputMultiplex, compound_map),
                       col_names = c("cpd", "formula", "mz", "rt", "rtwin", "batch_"))
  map_add <- map_data %>%
    filter(batch_ == batch) %>%
    rename(batch = batch_) %>%
    mutate(raw_file = raw_file, compound_map = compound_map)
  map_add

})

# multiplicate raw files

pwalk(tables, function(...) {
  data <- list(...)
  target_name <- data$cpd %>% str_replace_all(" ", "_")
  fs::file_copy(
    fs::path(inputMultiplex, data$raw_file),
    fs::path(inputDir, glue("{target_name}.mzML"))
  )
})

# Create compound list with help of Cosomi

cpdList <- pmap_dfr(tables, function(...) {
  data <- list(...)

  cpd_id <- data$cpd %>% str_match("([0-9]+)_.*") %>% `[`(,2)


  return(list(
    `UCHEM ID`=as.numeric(cpd_id),
    RT = data$rt
  ))

}) %>% left_join(cosomi_data)

cpdListRmb <- cpdList %>%
  transmute(
    ID = `UCHEM ID`,
    RT = RT,
    Compound = Name,
    SMILES=`Smiles str`,
    InChIKey = `Inchi key str`
  )

write_csv(
  cpdListRmb,
  compoundsCsv
)
