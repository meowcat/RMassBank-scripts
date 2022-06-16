
# Peak EIC correlation while filtering out zero rows
.eicScoreCor <- function(prec, eic) {
  eicSum <- colSums(abs(eic))
  .eicScore <- cor(prec, eic[,eicSum > 0,drop=FALSE])
  eicScore <- rep(NA_real_, ncol(eic))
  eicScore[eicSum > 0] <- .eicScore
  eicScore
} 

# Peak EIC dot product score
.eicScoreDot <- function(prec, eic) {
  normX <- sqrt(sum(prec^2))
  normY <- sqrt(colSums(eic^2))
  dot <- t(prec) %*% eic
  score <- dot / (normX * normY)
  score[normY == 0] <- NA
  score
}


f1 <- function(sens, spec)
  2 * sens * spec / (sens + spec)
fBeta <- function(beta) function(sens, spec)
  (1 + beta^2) * (sens * spec) / ((spec * beta^2) + sens)

checkDuplicates <- function(w) {
  id_adducts <- map_dfr(w@spectra %>% as.list(), function(cpd) {
    list(id = cpd@id, mode = cpd@mode, name = cpd@name, 
         acq_id = ifelse(
           "acq_ID" %in% names(attributes(cpd)),
           attr(cpd@id, "acq_ID"), 
           NA_character_)
    )
  })
  id_adducts %>% group_by(id, mode) %>% add_tally() %>% filter(n > 1) 
}


# cmap: adapt purrr::map at level 1. 
# c stands for compound, as opposed to smap, which works on s-pectra;
# note that cmap works on any SimpleList (also RmbSpectrum2List) and smap is like map_depth(2)
cmap <- function(x, ...) {
  x %>% as.list() %>% map(...) %>% as("SimpleList")
}
cmap_chr <- function(x, ...) {
  x %>% as.list() %>% map_chr(...)
}
cmap_int <- function(x, ...) {
  x %>% as.list() %>% map_int(...)
}
cmap_dbl <- function(x, ...) {
  x %>% as.list() %>% map_dbl(...)
}
cmap_lgl <- function(x, ...) {
  x %>% as.list() %>% map_lgl(...)
}
cmap_dfr <- function(x, ...) {
  x %>% as.list() %>% map_dfr(...)
}

smap <- function(x, ...) {
  x %>% as.list() %>% map(function(xx) {
    xx@children <- xx@children %>% cmap(...)
    xx
  }) %>% as("SimpleList")
}



