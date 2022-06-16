library(here)

library(tidyverse)
library(glue)


source(here("viewer.R"))

charge_str <- "pH"
backupPath <- glue("results/viewer_status_{charge_str}.RData")


w <- loadMsmsWorkspace(
  glue("results/spectra-{charge_str}-eics-score.RData"))

# Apply the specified specOK cutoff to reduce excessive work for nothing
w@spectra <- w@spectra %>% smap(function(ch) {
  d <- getData(ch)
  mx <- max(d$intensity, 0, na.rm = TRUE)
  ch@ok <- ch@ok & (mx > as.numeric(getOption("RMassBank")$filterSettings$specOkLimit))
  ch
})

# Check to exclude 0-spectra-found cpds
w@spectra <- w@spectra %>% cmap(function(cpd) {
  n_ok <- sum(cmap_lgl(cpd@children, ~ .x@ok))
  cpd@found <- n_ok > 0
  cpd
})

archiveResults(
  w, 
  glue("results/spectra-{charge_str}-autoreview.RData")
)


review <- viewer(w, backupPath = backupPath)

# Export review results as CSV
write_csv(tibble(
  cpd = seq_along(review$cpdOk),
  name = names(w@spectra),
  ok = review$cpdOk),
  file = glue("results/review_{charge_str}_cpd_ok.csv"))

spec_ok <-  pmap_dfr(
  tibble(.x = review$specOk %>% as.list(),
         .y = seq_along(review$specOk),
         ce = cmap(w@spectra, function(cpd) cmap_chr(cpd@children, ~.x@collisionEnergy)) %>% as.list(),
         name = names(w@spectra)), 
  function(.x, .y, name, ce) {
    tibble(cpd = .y, name = name, spectrum = seq_along(.x), ce = ce, ok = .x)
  })

write_csv(spec_ok, 
          file = glue("results/review_{charge_str}_spec_ok.csv"))

write_file(as.character(review$score_cutoff), 
           file = glue("results/review_{charge_str}_score_cutoff.csv"))
