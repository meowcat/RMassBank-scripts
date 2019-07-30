library(shiny)
library(plyr)
library(RMassBank)
library(MSnbase)

viewer <- function(w)
{
  ui <- fixedPage(
    sidebarLayout(
      sidebarPanel(
          selectInput("compound", "Compound:" ,names(w@spectra), size=20, selectize=FALSE)
      ),
      mainPanel(
        # first row: select cpd and MS

        tabsetPanel(
          tabPanel("spectra",
                   fixedRow(
                     numericInput("spectrum", "Spectrum:", 1,1,1)
                   ),
                   # second row: display MS
                   fixedRow(
                     column(7, plotOutput("spectrumPlot")),
                     column(5, tabsetPanel(
                       tabPanel("boxplot", plotOutput("eicCorPlot")),
                       tabPanel("xyplot", plotOutput("peaksXyPlot"))
                                 )
                   )),
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
      output$eicCorPlot <- renderPlot({
        cor <- property(cpd@children[[spectrum]], "cor")
        if(is.null(cor))
          plot.new()
        else {
          d <- getData(cpd@children[[spectrum]])
          # remove "bad" versions of the peak if there is a "good" one:
          d <- d[order(!d$good, abs(d$dppm)),,drop=FALSE]
          d <- d[!duplicated(d$mz),,drop=FALSE]
          par(mar=c(3,2,1,1)+0.1)
          boxplot(cor ~ good, data = d)
          bp <- boxplot(cor ~ good, data=d, horizontal=F, outline=FALSE)
          points(jitter(match(d$good, bp$names)), d$cor)
          # , 
          #        cex=d$intensity/max((d$intensity)) * 4,
          #        col=colorRamp(c("blue", "red"))(
          #          pmin(10,abs(d$dppm), na.rm = TRUE) / 10 
          #        ),
          #        pch=19)
          #         
          
          eps <- 0.01
          # plot(hist(cor, breaks=40))
          # lines(density(cor, bw="SJ"))
          sdCor <- sd(d$cor, na.rm = TRUE)
          abline(h=cpd@children[[spectrum]]@info$singlePointCor, col="red")
          abline(h=bp$stats[5,1]+eps, col="blue", lty=2)
        }
      })
      output$peaksXyPlot <- renderPlot(xyplot(cpd, spectrum))
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
  par(mar=c(3,2,1,1)+0.1)
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

xyplot <- function(cpd, spectrum)  {
  cor <- property(cpd@children[[spectrum]], "cor")
  if(is.null(cor))
    plot.new()
  else {
    d <- getData(cpd@children[[spectrum]])
    # remove "bad" versions of the peak if there is a "good" one:
    
    d <- d[order(!d$good, abs(d$dppm)),,drop=FALSE]
    d <- d[!duplicated(d$mz),,drop=FALSE]
    
    colGood <- c("TRUE" = "green", "FALSE" = "red")
    pchPeak <- c("TRUE" = 21, "FALSE" = 23)
    par(mar=c(3,2,1,1)+0.1)
    plot(pmin(abs(d$dppm), 10, na.rm=TRUE), d$cor,
         pch=pchPeak[as.character(d$good)], cex=log10(d$intensity)+0.2,
         xlim=c(0,10))
    
    
    abline(h=cpd@children[[spectrum]]@info$singlePointCor, col="red")
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