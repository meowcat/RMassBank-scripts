library(shiny)
library(plyr)
library(RMassBank)
library(MSnbase)

viewer <- function(w)
{
  ui <- fixedPage(
        # first row: select cpd and MS
        fixedRow(
          selectInput("compound", "Compound:" ,names(w@spectra))
        ),
        tabsetPanel(
          tabPanel("spectra",
                   fixedRow(
                       numericInput("spectrum", "Spectrum:", 1,1,1)
                     ),
                     # second row: display MS
                   fixedRow(
                       plotOutput("spectrumPlot")
                     ),
                     # third row: display peak table
                   fixedRow(
                       dataTableOutput("peaksTable")
                     )
                   ),
          tabPanel("eic",
                   fixedRow(
                     plotOutput("eicPlot")
                   )
                   
        )

      )
  )
  server <- function(input, output, session)
  {
    observe({
      cpd <- as.character(input$compound)
      if(cpd %in% names(w@spectra))
      {
        updateNumericInput(session, "spectrum", value = 1, 
                           min = 1, max = length(w@spectra[[cpd]]@children))
        output$eicPlot <- renderPlot(plotEic(w@spectra[[cpd]]))
      }
    }
      )
    observe({
      cpd <- input$compound
      spectrum <- input$spectrum
      
      if(!(cpd %in% names(w@spectra)))
        return()
      cpd <- w@spectra[[cpd]]
      if(!(spectrum < length(cpd@children)))
        return()
      # get the complete peak table to find the total range
      mzRange <- range(unlist(
        lapply(cpd@children, function(sp) range(mz(sp)))
      ))
      peaks <- getData(cpd@children[[spectrum]])
      output$peaksTable <- renderDataTable(peaks)
      output$spectrumPlot <- renderPlot({
        plotSpectrum(peaks, mzRange)
        spectrumLegend(cpd, cpd@children[[spectrum]])
        })
    })
  }
  shinyApp(server=server, ui=ui)
}

par(mfrow=c(1,1))
plotSpectrum <- function(peaks, mzRange=NA)
{
  xlim <- range(peaks$mz, mzRange, na.rm = TRUE)
  maxint <- max(peaks$intensity, na.rm=TRUE)
  peaks$intrel <- peaks$intensity/maxint*100
  plot.new()
  plot.window(xlim=xlim+c(-5,5), ylim=c(-10,120))
  axis(1)
  axis(2)
  points(intrel ~ mz, data=peaks, type='h', lwd=1, col="black")
  points(intrel ~ mz, data=peaks[peaks$good,,drop=FALSE],
         type='h', lwd=2, col="red")
}

plotEic <- function(cpd)
{
    eic <- attr(cpd, "eic")
    plot(intensity ~ rt, eic, type='l')
    if(length(cpd@children) > 0)
    {
      msData <- laply(cpd@children, function(ch)
        c(rt=ch@rt[[1]], tic=ch@tic[[1]]))
      lines(tic ~ rt, msData, col="red", type='h')
    }
}

spectrumLegend <- function(cpd, sp)
{
  leg <- c("name" = cpd@name,
    "precursor" = paste0(sp@precursorMz, " (", cpd@mz, ")"),
    "int" = format(max(intensity(sp)),scientific = TRUE, digits=1),
    "ce" = collisionEnergy(sp)
  )
  leg <- paste(leg, sep=": ")
  legend(
    "topleft",
    bty="n",
    legend=leg
  )
       
}