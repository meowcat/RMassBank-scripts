# https://github.com/rstudio/DT/issues/93#issuecomment-111001538
library(shiny)
library(DT)

shinyInput <- function(FUN, len, id, value, ...) {
  inputs = character(len)
  if(length(value) == 1)
    value <- rep(value, len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), value = value[[i]], label = NULL, ...))
  }
  inputs
}

# obtain the values of inputs
shinyValue <- function(input, id, len) {
  unlist(lapply(seq_len(len), function(i) {
    value <- input[[paste0(id, i)]]
    if (is.null(value)) NA else value
  }))
}

.preDrawCallback <- JS('function() { Shiny.unbindAll(this.api().table().node()); }')
.drawCallback <- JS('function() { Shiny.bindAll(this.api().table().node()); } ')

# 
# 
# shinyApp(
#   ui = fluidPage(DT::dataTableOutput('x1'), verbatimTextOutput('x2')),
#   
#   server = function(input, output) {
#     # create a character vector of shiny inputs
#     
#     
#     # a sample data frame
#     res = data.frame(
#       v1 = shinyInput(numericInput, 100, 'v1_', value = 0),
#       v2 = shinyInput(checkboxInput, 100, 'v2_', value = TRUE),
#       v3 = rnorm(100),
#       v4 = sample(LETTERS, 100, TRUE),
#       stringsAsFactors = FALSE
#     )
#     
#     # render the table containing shiny inputs
#     output$x1 = DT::renderDataTable(
#       res, server = FALSE, escape = FALSE, selection = 'single', options = list(
#         preDrawCallback = .preDrawCallback,
#         drawCallback = .drawCallback)
#       )
#     
#     # print the values of inputs
#     output$x2 = renderPrint({
#       data.frame(v1 = shinyValue(input, 'v1_', 100), v2 = shinyValue(input, 'v2_', 100))
#     })
#   }
# )