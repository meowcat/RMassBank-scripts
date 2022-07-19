library(tidyverse)
library(shiny)
library(plyr)
library(RMassBank)
library(MSnbase)
library(shinydust)
library(rhandsontable)
library(glue)
library(fs)
library(DT)
# library(keys)
source("viewer-include.R")

metric_set <- "eicScoreCor"

debug <- FALSE
debug_message <- function(...) {
  if(debug)
    message(...)
}

charge_str_select <- function(charge_strs) {
  ui <- fluidPage(
    selectInput("selection", "Review spectra:", charge_strs),
    actionButton("ok", "OK")
  )
  server <- function(input, output) {
    observeEvent(input$ok, {
      stopApp(returnValue = input$selection)
    })
  }
  vwr = dialogViewer('RMassBank', width = 300, height = 50)
  return(runApp(shinyApp(server=server, ui=ui, options=list(width=50, height=50))))
  
}


viewer <- function(w, backupPath = fs::path(getwd(), "viewer_status.RData"))
{
  if(is.list(w)) {
    w_ <- reactiveValues(
      w = w$w,
      specOk = w$specOk,
      cpdOk = w$cpdOk,
      cpd = 1,
      spec = 1,
      score_cutoff = ifelse(!is.null(w$score_cutoff), w$score_cutoff, attr(w$w, "eicScoreFilter")[[metric_set]])
    )
  } else {
  w_ <- reactiveValues(w=w,
                       specOk = cmap(
                         w@spectra, function(cpd)
                           cmap_lgl(cpd@children, function(sp)
                             sp@ok)),
                       cpdOk = unlist(lapply(
                         w@spectra, function(cpd)
                           cpd@found
                       )),
                       cpd = 1,
                       spec = 1,
                       score_cutoff = attr(w, "eicScoreFilter")[[metric_set]]
  )}
  
  hotkeys <- c("+", "-", ".")
  
  ui <- fluidPage(
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
        width = 4,
        fixedRow(rHandsontableOutput("compound", height = "400px") ), #, "Compound:" ,names(w@spectra), size=20, selectize=FALSE)
        fixedRow(rHandsontableOutput("spectrum")),
        fixedRow(actionButton("quit", "Quit")),
        fixedRow(actionButton("reset", "Reset")),
        fixedRow(actionButton("backup", "Backup")),
        fixedRow(actionButton("restore", "Restore"))
      ),
      mainPanel(
        width = 8,
        # first row: select cpd and MS

        tabsetPanel(
          tabPanel("spectra",
                   # fixedRow(
                   #   numericInput("spectrum", "Spectrum:", 1,1,1)
                   # ),
                   # second row: display MS
                   fixedRow(
                     column(5, plotOutput("spectrumPlot")),
                     column(3, plotOutput("eicCorPlot", click = "eicCorPlot_click")),
                     column(4, plotOutput("eicPeakPlot"))
                   ),
                   # third row: display peak table
                   fixedRow(
                     DT::dataTableOutput("peaksTable")
                   )
          ),
          tabPanel("eic",
                   fixedRow(
                     plotOutput("eicPlot")
                   )
                   
          )
          
        )
      )
    ),
    # useKeys(),
    # keysInput("keys", hotkeys),
  )
  server <- function(input, output, session)
  {
    
    
    
    #output$ms2table <- 
    compound <- reactive({
      names(w_$w@spectra)[compoundIndex()]
    })
    
    compoundIndex <- reactive({
      #input$compound_select$select$r
      debug_message(w_$cpd)
      w_$cpd
    })
    
    observeEvent(input$compound_select, {
      r <- input$compound_select$select$r
      if((r != w_$cpd) & (length(r) == 1))
        w_$cpd <- r
    })
    
    observeEvent(input$eicCorPlot_click, {
      w_$score_cutoff <- input$eicCorPlot_click$y
    })
    
    # observeEvent(input$keys, {
    #   message(input$keys)
    #   actions <- list(
    #     "+" = function() w_$cpd <- min(w_$cpd + 1, length(w_$w@spectra)),
    #     "-" = function() w_$cpd <- max(w_$cpd - 1, 1),
    #     "." = function() w_$cpdOk[[w_$cpd]] <- !(w_$cpdOk[[w_$cpd]])
    #   )
    #   if(input$keys %in% names(actions))
    #     actions[[input$keys]]()
    # })
    
    
    output$compound <- renderRHandsontable({
      specNames <- names(w_$w@spectra)
      specModes <- cmap_chr(w_$w@spectra, ~.x@mode)
      adductTable <- RMassBank:::getAdductInformation("")
      specAdducts <- adductTable[match(specModes, adductTable$mode), "adductString"]
      countOK <- cmap_int(w_$w@spectra, function(cpd) {
        sum(cmap_lgl(cpd@children, ~ .x@ok))
      })
      debug_message(glue("{length(specNames)} cpds, {length(w_$cpdOk)} list of cpds"))
      # message(length(specNames))
      countTot <- cmap_int(w_$w@spectra, ~ length(.x@children))
      df <- tibble(
        ok = w_$cpdOk,
        #ok = rep(TRUE, length(specNames)),
        name = specNames,
        adduct = specAdducts,
        num_ok = as.character(glue("{countOK}/{countTot}")))
      rhandsontable(df, 
                    selectCallback = TRUE, 
                    rowHeaders = NULL,
                    rowhighlight = 7) %>%
        hot_rows(rowHeights = 15) %>%
        # hot_cols(renderer = "
        # function(instance, td, row, col, prop, value, cellProperties) {
        #       Handsontable.renderers.TextRenderer.apply(this, arguments);
        #       tbl = this.HTMLWidgets.widgets[0];
        # 
        #       hrows = tbl.params.rowhighlight;
        #       hrows = hrows instanceof Array ? hrows : [hrows];
        #       if (hrows.includes(row)) {
        #         td.style.background = 'pink';
        #       }
        #       return td;
        # }
        #          ") %>%
        hot_col(col = "ok", type = "checkbox") %>%
        hot_col(col = "name", readOnly = TRUE) %>%
        hot_col(col = "adduct", readOnly = TRUE) %>%
        hot_col(col = "num_ok", readOnly = TRUE) %>%
        hot_table(highlightRow = TRUE)
    })
    
    
    # output$compound <-  DT::renderDataTable({
    #   specNames <- names(w_$w@spectra)
    #   specModes <- map_chr(w_$w@spectra %>% as.list(), ~.x@mode)
    #   adductTable <- RMassBank:::getAdductInformation("")
    #   specAdducts <- adductTable[match(specModes, adductTable$mode), "adductString"]
    #   data.frame(ok = shinyInput(checkboxInput, length(w_$w@spectra), "cpd_selected", 
    #                              value=isolate({w_$cpdOk})),
    #              name = glue("{specNames} {specAdducts}")) %>%
    #     datatable(selection = "single", filter="none",
    #               class="compact",
    #               options = list("dom" = "t", scrollY=350, paging=FALSE,
    #                              "preDrawCallback" = .preDrawCallback,
    #                              "drawCallback" = .drawCallback),
    #               escape = FALSE,
    #               rownames = FALSE) %>%
    #     DT::formatStyle(columns = c(1,2), fontSize = '75%')
    #     
    # })
    
    spectrum <- reactive({
      sel <- input$spectrum_select
      return(sel$select$r)
    })
    
    output$spectrum <- renderRHandsontable({
      cpd <- compoundIndex()
      if(length(cpd) == 0)
        return()
      if(between(cpd, 1, length(w_$w@spectra)))
      {
        isolate({
        debug_message(str(cpd))
        debug_message(glue("{w_$specOk[[cpd]]} "))
        df <- data.frame(
          ok = isolate({w_$specOk[[cpd]]}),
          id = seq_along(w_$w@spectra[[cpd]]@children),
          int = w_$w@spectra[[cpd]]@children %>% as.list() %>% map_dbl(~ getData(.x) %>% pull(intensity) %>% max()) %>% format(scientific = T, digits = 2)
        )
        rhandsontable(df, selectCallback = TRUE) %>%
          hot_col(col = "ok", type = "checkbox") %>%
          hot_col(col = "id", readOnly = TRUE) %>%
          hot_col(col = "int", readOnly = TRUE) %>%
          hot_table(highlightRow = TRUE)
        })
      }
        
      #   
      #   %>%
      #     datatable(selection = "single", filter="none",
      #               class="compact", extensions = "KeyTable",
      #               options = list(
      #                 "dom" = "tp", keys = TRUE,
      #                 "preDrawCallback" = .preDrawCallback,
      #                 "drawCallback" = .drawCallback),
      #               escape = FALSE,
      #               rownames = FALSE) %>%
      #     DT::formatStyle(columns = c(1,2), fontSize = '75%'))
      # }
    })
    
    
    
    observe({
      #
      df <- hot_to_r(input$spectrum)
      isolate({
        idx <- compoundIndex()  
        debug_message(glue("setting {idx}, old: {w_$specOk[[idx]]}, new: {df$ok} \n"))
        if(is.null(df$ok))
          debug_message("{w_$specOk[[idx]]} would be NULLed - we skip this")
        else(
          if(length(idx) > 0)
            w_$specOk[[idx]] <- df$ok  
        )
      })
      
    })

    observe({
      #out <<- input$compound
      df <- hot_to_r(input$compound)
      isolate({
        # cpd <- compound()
        # idx <- compoundIndex()
        # if(length(cpd) == 0)
        #   return()
        # if(length(idx) == 0)
        #   return()
        if(length(df$ok) == length(w_$cpdOk))
          w_$cpdOk <- df$ok
        else
          debug_message("incorrect length of cpdOK - not setting (yet?)")
      })
    })

    
    observe({
      cpd <- compound()
      output$eicPlot <- renderPlot(plotEic(w_$w@spectra[[cpd]]))
    })
    
    
    peaksData <- reactive({
      cpdName <- compound()
      spectrum <- spectrum()
      if(length(cpdName) == 0)
        return()
      if(length(spectrum) == 0)
        return()
      
      if(!(cpdName %in% names(w_$w@spectra)))
        return()
      cpd <- w_$w@spectra[[cpdName]]
      if(!(spectrum <= length(cpd@children)))
        return()
      # get the complete peak table to find the total range
      mzRange <- range(unlist(
        lapply(cpd@children, function(sp) range(mz(sp)))
      ))
      
      sp <- cpd@children[[spectrum]]
      peaks <- getData(sp)
      
      # Keep best entry per peak
      peaksFiltered <- peaks %>%
        rowid_to_column("rowid") %>%
        group_by(mz) %>%
        arrange(!good, desc(formulaMultiplicity), abs(dppm)) %>%
        slice(1)
      
      return(list(
        cpd=cpd,
        cpdName=cpdName,
        sp=sp,
        spectrum=spectrum,
        peaks=peaks,
        peaksFiltered=peaksFiltered,
        mzRange=mzRange
      ))
    })
    
    observe({
      
      data <- peaksData()
      if(length(data) == 0)
        return()
      
      peaks <- data$peaksFiltered
      cpd <- data$cpd
      spectrum <- data$spectrum
      mzRange <- data$mzRange
      cpdName <- data$cpdName
      sp <- data$sp
      
      peaks$mzPrecursor <- rep(cpd@mz, nrow(peaks))
      if(!("filterOK" %in% colnames(peaks)))
        peaks$filterOK <- rep(FALSE, nrow(peaks))
      if(!("formulaSource" %in% colnames(peaks)))
        peaks$formulaSource <- rep("", nrow(peaks))
      if(!("bestMultiplicity" %in% colnames(peaks)))
          peaks$bestMultiplicity <- rep(0, nrow(peaks))

      output$peaksTable <- renderDataTable({
        peaks %>%
          mutate(flag = glue(
            '{if_else(noise, "noise", "")}{if_else(low, "low ","")}{if_else(satellite, "sat ", "")}{if_else(good, "good ","")}{if_else(filterOK, "filter+", "")}')) %>%
          select(-low, -satellite, -good, -rawOK, -formulaSource, -dppmBest, -bestMultiplicity, -filterOK, -noise, -rowid) %>%
          dplyr::rename(cor = eicScoreCor, dot = eicScoreDot, int = intensity, `#mf` = formulaCount, `#spec` = formulaMultiplicity) %>%
          dplyr::relocate(flag, .before = formula) %>%
          dplyr::relocate(cor, .after = dppm) %>%
          dplyr::relocate(mzRaw, .before = mzPrecursor) %>%
          datatable(selection = "single") %>%
          formatRound(c("mz", "mzRaw", "mzCalc", "mzPrecursor"), digits = 4) %>%
          formatRound(c("dppm", "cor", "dot"), digits = 2) %>%
          formatSignif("int", 2) })
      output$spectrumPlot <- renderPlot({
        plotSpectrum(peaks, mzRange=mzRange, score_cutoff = w_$score_cutoff, metric = metric_set)
        spectrumLegend(cpd, sp, cpdName)
        })
      output$eicCorPlot <- renderPlot(plotEicCor(w_$score_cutoff, cpd, spectrum, metric = metric_set))
      #output$peaksXyPlot <- renderPlot(xyplot(cpd, spectrum))
    })
    
    
    eicData <- reactive({
      
      data <- peaksData()
      if(length(data) == 0)
        return()
      
      cpd <- data$cpd
      sp <- data$sp
      valid_rowids <- as.character(data$peaksFiltered$rowid)
      
      eicParent <- attr(cpd, "eic") %>% mutate(
        rowid = "parent",
        precursorScan = scan)
      eicPeak <- attr(sp, "eic") %>%
        bind_rows(.id = "rowid") %>%
        filter(rowid %in% valid_rowids)
      
      validPrecursorScans <- unique(eicPeak$precursorScan)
      
      return(bind_rows(list(
        eicParent %>% filter(precursorScan %in% validPrecursorScans),
        eicPeak
      )))
      
    })
    
    # Plot EIC for parent and selected peak
    observe({
      
      data <- peaksData()
      if(length(data) == 0)
        return()
      peakSelected <- input$peaksTable_rows_selected
      
      peaks <- data$peaksFiltered
      cpd <- data$cpd
      spectrum <- data$spectrum
      mzRange <- data$mzRange
      cpdName <- data$cpdName
      sp <- data$sp
      
      origRowId <- as.character(peaks$rowid[peakSelected])
      
      eic <- eicData() %>%
        group_by(rowid) %>%
        dplyr::mutate(relint = intensity / max(intensity, na.rm = TRUE))
        
      output$eicPeakPlot <- renderPlot({
        eic %>%
          filter(rowid %in% c("parent", origRowId)) %>%
          ggplot() +
          aes(x=rt, y=relint, color = rowid) +
          geom_line() + 
          geom_point() +
          theme_minimal()
      })
      
      message(peakSelected)
      message(origRowId)
      message(data$peaks$mz[origRowId])
      # 
    })
    
    observeEvent(input$quit, {
      stopApp(returnValue = reactiveValuesToList(w_))
    })
    
    observeEvent(input$reset, {
      debug_message("reset")
      if(is.null(input$reset))
        return()
      debug_message("reset executed")
      w_$specOk <- lapply(w_$w@spectra, function(cpd)
                               unlist(lapply(cpd@children, function(sp)
                                 sp@ok)))
      w_$cpdOk <- unlist(lapply(
                             w_$w@spectra, function(cpd)
                               cpd@found
                           ))
      })
    
    
    observeEvent(input$backup, {
      if(is.null(input$backup))
        return()
      e <- new.env()
      e$specOk <- w_$specOk
      e$cpdOk <- w_$cpdOk
      e$score_cutoff <- w_$score_cutoff
      save(list = c("specOk", "cpdOk", "score_cutoff"), envir = e, file=backupPath)
      showNotification("Backup complete")
    })
    
    
    
    observeEvent(input$restore, {
      if(is.null(input$restore))
        return()
      e <- new.env()
      load(backupPath, envir = e)
      w_$specOk <- e$specOk
      w_$cpdOk <- e$cpdOk
      if(!is.null(e$score_cutoff))
        w_$score_cutoff <- e$score_cutoff
      showNotification("Restore complete")
    })
    
  }
  return(runApp(shinyApp(server=server, ui=ui)))
}

par(mfrow=c(1,1))
plotSpectrum <- function(peaks, mzRange=NA, score_cutoff, metric = metric_set)
{
  xlim <- range(peaks$mz, mzRange, na.rm = TRUE)
  maxint <- max(peaks$intensity, na.rm=TRUE)
  # select only the best peak
  peaks <- peaks %>% group_by(mz) %>% arrange(!good, abs(dppm)) %>% slice(1)
  peaks$intrel <- peaks$intensity/maxint*100
  peaks$corOk <- peaks[,metric] > score_cutoff
  par(mar=c(3,2,1,1)+0.1)
  plot.new()
  plot.window(xlim=xlim+c(-5,5), ylim=c(-10,120))
  axis(1)
  axis(2)
  points(intrel ~ mz, data=peaks, type='h', lwd=1, col="black")

  points(intrel ~ mz, data=peaks[peaks$good,,drop=FALSE],
         type='h', lwd=2, col="red")
  points(intrel ~ mz, data=peaks[peaks$corOk & peaks$good,,drop=FALSE],
         type='h', lwd=2, col="blue")
  points(intrel ~ mz, data=peaks[peaks$corOk & !peaks$good,,drop=FALSE],
         type='h', lwd=2, col="green")
}

plotEicCor <- function(score_cutoff, cpd, spectrum, metric = metric_set) {
  cor <- property(cpd@children[[spectrum]], metric)
  if(is.null(cor))
    plot.new()
  else {
    d <- getData(cpd@children[[spectrum]])
    
    # Keep best entry per peak
    d <- d %>%
      dplyr::group_by(mz) %>%
      dplyr::arrange(!good, desc(formulaMultiplicity), abs(dppm)) %>%
      dplyr::slice(1)
    # 
    # # remove "bad" versions of the peak if there is a "good" one:
    # d <- d[order(!d$good, abs(d$dppm)),,drop=FALSE]
    # d <- d[!duplicated(d$mz),,drop=FALSE]
    par(mar=c(3,2,1,1)+0.1)
    d <- d %>% 
      dplyr::mutate(metric = replace_na(!!sym(col), 0))
    d$good <- if_else(d$good, "MF found", "no MF") %>% factor(levels = c("no MF", "MF found"))
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
    
    abline(h=score_cutoff, col="red")
    title(ylab = "correlation")
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



spectrumLegend <- function(cpd, sp, cpdName)
{
  leg <- c("name" = glue("{cpd@name} ({cpdName})"),
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