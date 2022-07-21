library(tidyverse)
library(shiny)
library(shinyWidgets)
library(plyr)
library(RMassBank)
library(MSnbase)
library(shinydust)
library(rhandsontable)
library(glue)
library(fs)
library(DT)
library(zoo)
# library(keys)
source("viewer-include.R")

metric_set <- "eicScoreCor"

colors_threshold <- list(
  "global" = "red",
  "cpd" = "blue",
  "spec" = "grey"
)
html_threshold = Vectorize(function(option) {
  HTML(glue("<span style='color: {colors_threshold[option]}'> {option} </span>"))
})

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

# https://laustep.github.io/stlahblog/posts/DTcallbacks.html#select-rows-on-click-and-drag
# Note: "key-focus" is the event triggered because KeyTable key nav *already* causes to
# move on the grid, and triggers focusing of cells.
# In contrast, "key" would only handle events that are *not* handled by key-focus or other
# built-in KeyTable functionality.
js_select_dt <- c(
  "var dt = table.table().node();",
  "var tblID = $(dt).closest('.datatables').attr('id');",
  "var inputName = tblID + '_rows_selected'",
  "var incrementName = tblID + '_rows_selected2_increment'",
  "table.on('key-focus', function(e, datatable, cell, originalEvent){",
  "   console.log(originalEvent);",
  #"   console.log(key);",
  "  if (originalEvent.type === 'keydown'){",
  #"  if (keysReact.indexOf(key) > -1){",
  "    table.rows().deselect(); ",
  "    table.row(cell[0][0].row).select();",
  "    row = table.rows({selected: true})",
  #"    Shiny.setInputValue(inputName, [row[0]]);",
  # Note: this ID is zero-based so add one
  "    Shiny.setInputValue(inputName, [parseInt(row[0]) + 1]);",
  "  }",
  "});"
)

generateSpecOk <- function(w) {
  cmap(
    w@spectra, function(cpd)
      tibble(
        ok = cmap_lgl(cpd@children, function(sp)
          sp@ok),
        threshold = NA_real_
      )
  )
}

generateCpdOk <- function(w) {
  tibble(
    ok = unlist(lapply(
      w@spectra, function(cpd)
        cpd@found
    )),
    threshold = NA_real_)
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
                       specOk = generateSpecOk(w),
                       cpdOk = generateCpdOk(w),
                       cpd = 1,
                       spec = 1,
                       score_cutoff = attr(w, "eicScoreFilter")[[metric_set]]
  )}
  
  frozen <- reactiveVal(FALSE)
  
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
                     column(2,
                            radioGroupButtons(
                              inputId = "thresholdMode",
                              label = "Threshold",
                              choiceNames = html_threshold(names(colors_threshold)) %>% set_names(NULL),
                              choiceValues = names(colors_threshold) %>% set_names(NULL),
                              size = "xs"
                            ),
                            actionLink("clearThreshold", "clear"),
                            plotOutput("eicCorPlot", click = "eicCorPlot_click", height = 350),
                            #tags$style("#thresholdMode {font-size: 10px;}"),

                            ),
                     column(5, plotOutput("eicPeakPlot"))
                   ),
                   fixedRow(
                     prettyCheckbox("filterPeaksTable", "Show only filtered peaks")
                   ),
                   # third row: display peak table
                   fixedRow(
                     DT::DTOutput("peaksTable")
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
      score_cutoff_ <- input$eicCorPlot_click$y
      isolate({
        idx <- compoundIndex()  
        spIdx <- spectrum()
        # set the threshold either per cpd, spectrum or globally
        switch (input$thresholdMode,
          "global" = {
            w_$score_cutoff <- score_cutoff_
          },
          "cpd" = {
            cutoffs <- w_$cpdOk$threshold
            cutoffs[[idx]] <- score_cutoff_
            w_$cpdOk$threshold <- cutoffs
          },
          "spec" = {
            # message(glue("th{score_cutoff_}"))
            cutoffs <- w_$specOk[[idx]]$threshold
            cutoffs[[spIdx]] <- score_cutoff_
            # message(cutoffs)
            w_$specOk[[idx]]$threshold <- cutoffs
          }
        ) # switch
      }) # isolate
    }) # observeEvent
    
    observeEvent(input$clearThreshold, {
      score_cutoff_ <- NA_real_
      isolate({
        idx <- compoundIndex()  
        spIdx <- spectrum()
        # reset the threshold either per cpd, spectrum or globally
        switch (input$thresholdMode,
                "global" = {
                  w_$score_cutoff <- attr(w_$w, "eicScoreFilter")[[metric_set]]
                },
                "cpd" = {
                  cutoffs <- w_$cpdOk$threshold
                  cutoffs[[idx]] <- score_cutoff_
                  w_$cpdOk$threshold <- cutoffs
                },
                "spec" = {
                  # message(glue("th{score_cutoff_}"))
                  cutoffs <- w_$specOk[[idx]]$threshold
                  cutoffs[[spIdx]] <- score_cutoff_
                  # message(cutoffs)
                  w_$specOk[[idx]]$threshold <- cutoffs
                }
        ) # switch
      }) # isolate
    }) # observeEvent
    
    
    # collect currently set thresholds and the active one for the current
    # spectrum and cpd. Specific overrules general.
    validThreshold <- reactive({
      idx <- compoundIndex()  
      spIdx <- spectrum()
      th <- list(
        global = w_$score_cutoff,
        cpd = w_$cpdOk$threshold[[idx]],
        spec = w_$specOk[[idx]]$threshold[[spIdx]]
      )
      which_valid <- tail(which(!is.na(th)), 1)
      th[["which_valid"]] <- names(th)[[which_valid]]
      th[["valid"]] <- th[[which_valid]]
      
      th
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
    
    
    output$compound <- renderRHandsontable({ withProgress( message="wait", {
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
        ok = w_$cpdOk$ok,
        #ok = rep(TRUE, length(specNames)),
        name = specNames,
        adduct = specAdducts,
        threshold = w_$cpdOk$threshold,
        num_ok = as.character(glue("{countOK}/{countTot}")))
      rh <- rhandsontable(df, 
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
        hot_col(col = "threshold", readOnly = TRUE) %>%
        hot_table(highlightRow = TRUE)
      frozen(FALSE)
      rh
    }) # withProgress
      }) # renderRHandsonTable
    
    
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
      
      
      ok <- w_$specOk
      if(between(cpd, 1, length(w_$w@spectra)))
      {
        isolate({
        # observe w_$w
        
        
          
        debug_message(str(cpd))
        debug_message(glue("{w_$specOk[[cpd]]} "))
        df <- data.frame(
          ok = isolate({w_$specOk[[cpd]]$ok}),
          id = seq_along(w_$w@spectra[[cpd]]@children),
          int = w_$w@spectra[[cpd]]@children %>% as.list() %>% map_dbl(~ getData(.x) %>% pull(intensity) %>% max()) %>% format(scientific = T, digits = 2),
          threshold = isolate({w_$specOk[[cpd]]$threshold})
        )
        rhandsontable(df, selectCallback = TRUE) %>%
          hot_col(col = "ok", type = "checkbox") %>%
          hot_col(col = "id", readOnly = TRUE) %>%
          hot_col(col = "int", readOnly = TRUE) %>%
          hot_col(col = "threshold", readOnly = TRUE) %>%
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
    
    
    # Update the spectra checked status (the "ok" column in the dataframe
    # retrieved from HandsOnTable, which is the checkboxes) into the w_$specOK
    # variable, which will finally be exported as a textfile and decide which
    # spectra are exported
    observe({
      #
      #out <<- input$compound
      if(frozen())
        return()
      df <- hot_to_r(input$spectrum)
      isolate({
        idx <- compoundIndex()  
        debug_message(glue("setting {idx}, old: {w_$specOk[[idx]]$ok}, new: {df$ok} \n"))
        if(is.null(df$ok))
          debug_message("{w_$specOk[[idx]]} would be NULLed - we skip this")
        else(
          if(length(idx) > 0)
            w_$specOk[[idx]]$ok <- df$ok  
        )
      })
      
    })

    
    # Update the compound checked status (the "ok" column in the dataframe
    # retrieved from HandsOnTable, which is the checkboxes) into the w_$cpdOk
    # variable, which will finally be exported as a textfile and decide which
    # compounds are exported
    observe({
      #out <<- input$compound
      if(frozen())
        return()
      frozen(TRUE)
      df <- hot_to_r(input$compound)
      isolate({
        # cpd <- compound()
        # idx <- compoundIndex()
        # if(length(cpd) == 0)
        #   return()
        # if(length(idx) == 0)
        #   return()
        if(length(df$ok) == length(w_$cpdOk$ok))
          w_$cpdOk$ok <- df$ok
        else
          debug_message("incorrect length of cpdOK - not setting (yet?)")
      })
    })

    # Plot EIC on second tab. This is somewhat obsolete with the new combined
    # EIC plot including fragments.
    observe({
      cpd <- compound()
      output$eicPlot <- renderPlot(plotEic(w_$w@spectra[[cpd]]))
    })
    
    
    # Collect all data relating to the selected spectrum (cpd + child):
    # the RMassBank objects, peak table, deduplicated peak table, extras
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
    
    # Render peak table and correlations boxplot
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

      validTh <- validThreshold()
      filtered <- input$filterPeaksTable
      
      output$peaksTable <- DT::renderDT({
        
        
        peaks %>%
          dplyr::mutate(filterOK = !!sym(metric_set) >= validTh$valid) %>%
          dplyr::filter(!filtered | filterOK) %>%
          dplyr::mutate(flag = glue(
            '{if_else(noise, "noise", "")}{if_else(low, "low ","")}{if_else(satellite, "sat ", "")}{if_else(good, "good ","")}{if_else(filterOK, "filter+", "")}')) %>%
          dplyr::select(-low, -satellite, -good, -rawOK, -formulaSource, -dppmBest, -bestMultiplicity, -filterOK, -noise, -rowid) %>%
          dplyr::rename(cor = eicScoreCor, dot = eicScoreDot, int = intensity, `#mf` = formulaCount, `#spec` = formulaMultiplicity) %>%
          dplyr::relocate(flag, .before = formula) %>%
          dplyr::relocate(cor, .after = dppm) %>%
          dplyr::relocate(mzRaw, .before = mzPrecursor) %>%
          datatable(
            # This datatable uses both shiny's select for conventional selection
            # and keytable + select for selection by keyboard (callback js_select_dt, see above).
            # The keyboard-selected row just overwrites the regular input$peaksTable_rows_selected
            # field.
            # Note: server=FALSE is required (see below) or this will return wrong row IDs
            # when using client-side sorting.
            selection = "single",
            editable = FALSE, 
            callback = JS(js_select_dt),
            extensions = c("KeyTable", "Select"),
            options = list(
              keys = TRUE,
              select = TRUE
            )
            ) %>%
          formatRound(c("mz", "mzRaw", "mzCalc", "mzPrecursor"), digits = 4) %>%
          formatRound(c("dppm", "cor", "dot"), digits = 2) %>%
          formatSignif("int", 2) }, 
        server = FALSE) # Note: server = FALSE is required (see above)
      output$spectrumPlot <- renderPlot({
        plotSpectrum(peaks, mzRange=mzRange, score_cutoff = validTh$valid, metric = metric_set)
        spectrumLegend(cpd, sp, cpdName)
        })
      output$eicCorPlot <- renderPlot(plotEicCor(validTh, cpd, spectrum, metric = metric_set))
      #output$eicCorPlot <- renderPlot(plotEicCor(1, cpd, spectrum, metric = metric_set))
      #output$peaksXyPlot <- renderPlot(xyplot(cpd, spectrum))
    })
    
    
    # Collect EIC from parent and fragments into longform tibble
    # from the attr(cpd, "eic") (a single chromatogram) and the 
    # attr(sp, "eic") (a list of chromatograms, one for each peak of getData())
    # Then deduplicate but keep the original rowid to know which peak in 
    # peaksFiltered belongs to the eic.
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
      
      # get selected row
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
        filter(rowid %in% c("parent", origRowId)) %>%
        dplyr::mutate(relint = intensity / max(intensity, na.rm = TRUE))
      

      
      getPeakLabel <- Vectorize(function(rowid_) {
        if(rowid_ == "parent")
          return(glue("parent\n(tic {format(sp@tic, digits=2, scientific=TRUE)})"))
        peaks %>% 
          filter(rowid == as.integer(rowid_)) %>%
          glue_data("{round(mz, 4)}\n({format(intensity, digits=2, scientific=TRUE)})")
        })
      
      # scaleInfo <- eic %>%
      #   dplyr::summarize(imax = max(intensity)) %>%
      #   glue_data("{getPeakLabel(rowid)}: int {format(imax, scientific=TRUE)}")
        
      output$eicPeakPlot <- renderPlot({
        eic %>%
          left_join(peaks %>% mutate(rowid = as.character(rowid)), by = "rowid") %>%
          ggplot() +
          aes(x=rt, y=relint, color = rowid) +
          scale_color_discrete(labels = getPeakLabel) +
          geom_line() + 
          geom_point() +
          # annotate("text", 
          #          label = scaleInfo %>% paste0(collapse="\n"),
          #          x = Inf, y = Inf,
          #          hjust = 1, vjust = 1) +
          theme_minimal()
        
      })
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
      w_$specOk <- generateSpecOk(w_$w)
      w_$cpdOk <- generateCpdOk(w_$w)
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

plotEicCor <- function(thresholds, cpd_, spectrum_, metric = metric_set) {
  cpd <- cpd_
  spectrum <- spectrum_
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
    metric_col <- metric
    d <- d %>% 
      dplyr::mutate(metric = replace_na(!!sym(metric_col), 0))
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
    
    for(threshold in names(colors_threshold)) {
      abline(h=thresholds[[threshold]], col=colors_threshold[[threshold]],
             lwd = if_else(threshold == thresholds[["which_valid"]], 2, 1))
    }
    
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
    "precursor" = glue("{round(sp@precursorMz,4)} (theor. {round(cpd@mz, 4)})"),
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