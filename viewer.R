library(shiny)
library(plyr)
library(RMassBank)
library(MSnbase)

source("viewer-include.R")

viewer <- function(w)
{
  w_ <- reactiveValues(w=w,
                       specOk = lapply(
                         w@spectra, function(cpd)
                           unlist(lapply(cpd@children, function(sp)
                             sp@ok))
                       ),
                       cpdOk = unlist(lapply(
                         w@spectra, function(cpd)
                           cpd@found
                       ))
                       )
  ui <- fixedPage(
    tags$style("
      .checkbox { /* checkbox is a div class*/
        line-height: 50%;
        margin-top: 2px;
        margin-bottom: 2px; /*set the margin, so boxes don't overlap*/
      }
      .form-group {
        margin-bottom: 0px;
      }
      "),
    sidebarLayout(
      sidebarPanel(
        fixedRow(dataTableOutput("compound")), #, "Compound:" ,names(w@spectra), size=20, selectize=FALSE)
        fixedRow(dataTableOutput("spectrum")),
        fixedRow(actionButton("quit", "Quit"))
      ),
      mainPanel(
        # first row: select cpd and MS

        tabsetPanel(
          tabPanel("spectra",
                   # fixedRow(
                   #   numericInput("spectrum", "Spectrum:", 1,1,1)
                   # ),
                   # second row: display MS
                   fixedRow(
                     column(7, plotOutput("spectrumPlot")),
                     column(5, plotOutput("eicCorPlot"))
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
    )
  )
  server <- function(input, output, session)
  {
    #output$ms2table <- 
    compound <- reactive({
      names(w_$w@spectra)[input$compound_rows_selected]
    })
    
    compoundIndex <- reactive({
      input$compound_rows_selected
    })
    
    output$compound <-  DT::renderDataTable({
      data.frame(ok = shinyInput(checkboxInput, length(w_$w@spectra), "cpd_selected", 
                                 value=isolate({w_$cpdOk})),
                 name = names(w_$w@spectra)) %>%
        datatable(selection = "single", filter="none",
                  class="compact",
                  options = list("dom" = "t", scrollY=350, paging=FALSE,
                                 "preDrawCallback" = .preDrawCallback,
                                 "drawCallback" = .drawCallback),
                  escape = FALSE,
                  rownames = FALSE) %>%
        DT::formatStyle(columns = c(1,2), fontSize = '75%')
        
    })
    
    spectrum <- reactive({
      input$spectrum_rows_selected
    })
    
    output$spectrum <-  DT::renderDataTable({
      cpd <- as.character(compound())
      if(length(cpd) == 0)
        return()
      if(cpd %in% names(w_$w@spectra))
      {
        return(data.frame(
          ok = shinyInput(checkboxInput, length(w_$w@spectra[[cpd]]@children),
                          paste0("spec_", compoundIndex(), "_selected"), 
                          isolate({w_$specOk[[compoundIndex()]]})),
          id = seq_along(w_$w@spectra[[cpd]]@children)
        ) %>%
                 datatable(selection = "single", filter="none",
                           class="compact", extensions = "KeyTable",
                           options = list(
                             "dom" = "tp", keys = TRUE,
                             "preDrawCallback" = .preDrawCallback,
                             "drawCallback" = .drawCallback),
                           escape = FALSE,
                           rownames = FALSE) %>%
          DT::formatStyle(columns = c(1,2), fontSize = '75%'))
      }
    })
    
    observe({
      cpd <- compound()
      spectrum <- spectrum()
      idx <- compoundIndex()
      if(length(cpd) == 0)
        return()
      if(length(spectrum) == 0)
        return()
      if(length(idx) == 0)
        return()
      
      nspec <- length(w_$w@spectra[[cpd]]@children)
      specSelected <- shinyValue(input,
                                 paste0("spec_", compoundIndex(), "_selected"),
                                 nspec)
      if(!any(is.na(specSelected))) {
        specOk <- w_$specOk
        specOk[[compoundIndex()]] <- specSelected
        w_$specOk  <- specOk
      }
      })
    
    observe({
      cpd <- compound()
      idx <- compoundIndex()
      if(length(cpd) == 0)
        return()
      if(length(idx) == 0)
        return()
      
      ncpd <- length(w_$w@spectra)
      cpdSelected <- shinyValue(input,
                                 paste0("cpd_selected"),
                                 ncpd)
      w_$cpdOk <- cpdSelected
    })
      
    
    observe({
      cpd <- compound()
      output$eicPlot <- renderPlot(plotEic(w_$w@spectra[[cpd]]))
    })
    
    observe({
      cpd <- compound()
      spectrum <- spectrum()
      if(length(cpd) == 0)
        return()
      if(length(spectrum) == 0)
        return()
      
      if(!(cpd %in% names(w@spectra)))
        return()
      cpd <- w_$w@spectra[[cpd]]
      if(!(spectrum <= length(cpd@children)))
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
      output$eicCorPlot <- renderPlot(plotEicCor(w_$w, cpd, spectrum))
      #output$peaksXyPlot <- renderPlot(xyplot(cpd, spectrum))
    })
    
    observeEvent(input$quit, {
      stopApp(returnValue = reactiveValuesToList(w_))
    })
  }
  return(runApp(shinyApp(server=server, ui=ui)))
}

par(mfrow=c(1,1))
plotSpectrum <- function(peaks, mzRange=NA)
{
  xlim <- range(peaks$mz, mzRange, na.rm = TRUE)
  maxint <- max(peaks$intensity, na.rm=TRUE)
  peaks$intrel <- peaks$intensity/maxint*100
  par(mar=c(3,2,1,1)+0.1)
  plot.new()
  plot.window(xlim=xlim+c(-5,5), ylim=c(-10,120))
  axis(1)
  axis(2)
  points(intrel ~ mz, data=peaks, type='h', lwd=1, col="black")
  points(intrel ~ mz, data=peaks[peaks$good,,drop=FALSE],
         type='h', lwd=2, col="red")
}

plotEicCor <- function(w, cpd, spectrum, metric = "eicScoreCor") {
  cor <- property(cpd@children[[spectrum]], metric)
  if(is.null(cor))
    plot.new()
  else {
    d <- getData(cpd@children[[spectrum]])
    # remove "bad" versions of the peak if there is a "good" one:
    d <- d[order(!d$good, abs(d$dppm)),,drop=FALSE]
    d <- d[!duplicated(d$mz),,drop=FALSE]
    par(mar=c(3,2,1,1)+0.1)
    d$metric <- d[,metric] %>% replace_na(0)
    #boxplot(metric ~ good, data = d)
    bp <- boxplot(metric ~ good, data=d, horizontal=F, outline=FALSE)
    points(jitter(match(d$good, bp$names)), d$metric)
    # , 
    #        cex=d$intensity/max((d$intensity)) * 4,
    #        col=colorRamp(c("blue", "red"))(
    #          pmin(10,abs(d$dppm), na.rm = TRUE) / 10 
    #        ),
    #        pch=19)
    #         
    
    abline(h=attr(w, "eicScoreFilter")[[metric]], col="red")
  }
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