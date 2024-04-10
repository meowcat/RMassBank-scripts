library(here)

library(tidyverse)
library(glue)
library(RMassBank)
loadRmbSettings("input/RmbSettings.ini")


source(here::here("viewer.R"))
source("functions.R")

charge_strs <- c("pH", "mH")
charge_str <- charge_str_select(charge_strs)
# charge_str <- "pH"
backupPath <- glue("results/viewer_status_{charge_str}.RData")


w <- loadMsmsWorkspace(
  glue("results/spectra-{charge_str}-eics-score.RData"))

# Remove completely empty cpds
not_empty <- w@spectra %>% cmap_lgl( ~ length(.x@children) > 0)

message(glue("{w@files[!not_empty]}: removing empty compound\n\n"))

w@spectra <- w@spectra[not_empty]
w@files <- w@files[not_empty]


# Apply the specified specOK cutoff to reduce excessive work for nothing
w@spectra <- w@spectra %>% smap(function(ch) {
  d <- getData(ch)
  mx <- max(d$intensity, 0, na.rm = TRUE)
  ch@ok <- ch@ok & (mx > as.numeric(getOption("RMassBank")$filterSettings$specOkLimit))
  ch
})


# Check to exclude 0-spectra-found cpds
w@spectra <- w@spectra %>% cmap(function(cpd) {
  n_ok <- sum(cmap_lgl(cpd@children, ~ isTRUE(.x@ok)))
  cpd@found <- n_ok > 0
  cpd
})


not_problematic <- w@spectra %>% cmap_lgl(function(cpd) isTRUE(cpd@found))
message(glue("{w@files[!not_problematic]}: removing problematic compound\n\n"))
w@spectra <- w@spectra[not_problematic]
w@files <- w@files[not_problematic]


archiveResults(
  w, 
  glue("results/spectra-{charge_str}-autoreview.RData")
)


review <- viewer(w, backupPath = backupPath)

# Export review results as CSV
write_csv(tibble(
  cpd = seq_along(review$cpdOk$ok),
  name = names(w@spectra),
  ok = review$cpdOk$ok,
  threshold = review$cpdOk$threshold),
  file = glue("results/review_{charge_str}_cpd_ok.csv"))

spec_ok <-
  tibble(.x = review$specOk %>% as.list(),
         .y = seq_along(review$specOk),
         ce = cmap(w@spectra, function(cpd) cmap_chr(cpd@children, ~.x@collisionEnergy)) %>% as.list(),
         name = names(w@spectra)) %>%
  pmap_dfr(
  function(.x, .y, name, ce) {
    tibble(cpd = .y, name = name, spectrum = seq_along(.x$ok), ce = ce, "ok" = .x$ok, "threshold" = .x$threshold)
  })

write_csv(spec_ok, 
          file = glue("results/review_{charge_str}_spec_ok.csv"))

write_file(as.character(review$score_cutoff), 
           file = glue("results/review_{charge_str}_score_cutoff.csv"))

