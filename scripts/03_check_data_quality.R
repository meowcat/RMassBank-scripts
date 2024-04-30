library(knitr)
library(rmarkdown)

knitr::opts_knit$set(
  root.dir = getwd()
)
rmarkdown::render(
  "include/check_data_quality.Rmd",
  output_dir="results"
)