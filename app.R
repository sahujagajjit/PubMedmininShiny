# Developer: JeeT
# A tool for PubMed data analysis for scientometrcs and text mining
# Date created: July 4, 2025
# Date last updated: July 4, 2025



# === Auto-install Required Packages ===
required_cran <- c(
  "shiny",
  "shinyFiles",
  "shinydashboard",
  "DT",
  "bibliometrix",
  "reshape2",
  "igraph",
  "visNetwork",
  "colourpicker",
  "shinycssloaders",
  "webshot2"
)

# Install CRAN packages if missing
install_if_missing <- function(pkgs) {
  new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
  invisible(lapply(pkgs, require, character.only = TRUE))
}

install_if_missing(required_cran)

# Check and install Chromium if not found (for webshot2 screenshots)
if (!webshot2::is_chromium_installed()) {
  webshot2::install_chromium()
}



ui <- dashboardPage(
  dashboardHeader(disable = TRUE),  # Disable default header
  dashboardSidebar(
    radioButtons("inputMode", "Input Mode:",
                 choices = c("Single File" = "file", "Directory" = "dir"),
                 selected = "file"),
    
    conditionalPanel(
      condition = "input.inputMode == 'file'",
      radioButtons("fileType", "File Type:",
                   choices = c("PubMed Format" = "pubmed", "CSV Format" = "csv"),
                   selected = "pubmed"),
      
      fileInput("singleFile", "Select File", accept = c(".txt", ".csv")),
      verbatimTextOutput("singleFilePath")
    ),
    
    conditionalPanel(
      condition = "input.inputMode == 'dir'",
      shinyDirButton("inputDir", "Select Input Directory", "Choose directory"),
      verbatimTextOutput("inputPath")
    ),
    
    actionButton("submit", "Submit")
  ),
  dashboardBody(
    # Custom header row with logos and title
    tags$head(
      tags$style(HTML("
        .custom-header {
          background-color: #3c8dbc;
          padding: 10px;
          color: white;
        }
        .custom-header img {
          height: 30px;
        }
        .custom-header-title {
          font-size: 24px;
          font-weight: bold;
          text-align: center;
          line-height: 50px;
        }
      "))
    ),
    fluidRow(
      class = "custom-header",
      column(2, tags$img(src = "GyanArrasLogo.png", align = "left")),  # Replace with your actual logo file
      column(8, div(class = "custom-header-title", "📚 PubMedMininShiny")),
      column(2, tags$img(src = "BN2_logo.jpg", align = "right"))  # Replace with your actual logo file
    ),
    uiOutput("mainUI"),
    div(class = "footer",
        tags$hr(),
        p("© ", format(Sys.Date(), "%Y"), " Bioinformatics Nexus Network (BN²), GyanArras Academy Pvt. Ltd. All rights reserved.",
          style = "text-align:center; color: #888; font-size: 14px; margin-bottom: 10px;")
    )
    
  )
)


server <- function(input, output, session) {
  volumes <- c(Home = fs::path_home(), "Root" = "/")
  
  inputDirPath <- reactiveVal(NULL)
  bibliodata <- reactiveVal(NULL)
  mergedPubmedTxtFile <- reactiveVal(NULL)
  rv <- reactiveValues(biblioResults = NULL, biblioSummary = NULL)
  
  # Directory selection
  shinyDirChoose(input, "inputDir", roots = volumes, session = session)
  observeEvent(input$inputDir, {
    inputDirPath(parseDirPath(volumes, input$inputDir))
  })
  
  output$inputPath <- renderText({
    req(inputDirPath())
    paste("Input Directory:", inputDirPath())
  })
  
  output$singleFilePath <- renderText({
    req(input$singleFile)
    paste("Selected File:", input$singleFile$name)
  })
  

  
  
  # Submit and process
  observeEvent(input$submit, {
    req(input$inputMode)
    
    withProgress(message = "Processing bibliometric data...", value = 0.3, {
      
      # Single File Mode
      if (input$inputMode == "file") {
        req(input$singleFile)
        filepath <- input$singleFile$datapath
        filename <- input$singleFile$name
        
        if (input$fileType == "pubmed" && !grepl("\\.txt$", filename, ignore.case = TRUE)) {
          showModal(modalDialog(
            title = "Invalid File Type",
            "Please upload a .txt file for PubMed format.",
            easyClose = TRUE,
            footer = modalButton("OK")
          ))
          return()
        }
        
        if (input$fileType == "csv" && !grepl("\\.csv$", filename, ignore.case = TRUE)) {
          showModal(modalDialog(
            title = "Invalid File Type",
            "Please upload a .csv file for CSV format.",
            easyClose = TRUE,
            footer = modalButton("OK")
          ))
          return()
        }
        
        if (input$fileType == "pubmed") {
          incProgress(0.5, detail = "Converting PubMed text to data frame...")
          df <- convert2df(file = filepath, dbsource = "pubmed", format = "plaintext")
        } else {
          incProgress(0.5, detail = "Reading CSV file...")
          df <- read.csv(filepath, stringsAsFactors = FALSE)
        }
        
        bibliodata(df)
      }
      
      # Directory Mode
      if (input$inputMode == "dir") {
        req(inputDirPath())
        incProgress(0.2, detail = "Reading and merging files...")
        files <- list.files(inputDirPath(), full.names = TRUE)
        
        if (length(files) > 1) {
          all_text <- sapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"))
          combined_text <- paste(all_text, collapse = "\n")
          tmpfile <- tempfile(fileext = ".txt")
          writeLines(combined_text, tmpfile)
          mergedPubmedTxtFile(tmpfile)
        } else {
          mergedPubmedTxtFile(NULL)
        }
        
        incProgress(0.4, detail = "Converting to data frame...")
        df <- convert2df(file = if (!is.null(mergedPubmedTxtFile())) mergedPubmedTxtFile() else files[1],
                         dbsource = "pubmed", format = "plaintext")
        bibliodata(df)
      }
      
      # Bibliometric analysis
      incProgress(0.1, detail = "Analyzing...")
      rv$biblioResults <- biblioAnalysis(bibliodata(), sep = ";")
      rv$biblioSummary <- summary(rv$biblioResults, k = 10, pause = FALSE)
    })
  })
  
  # Reactive: Year-wise data
  yearwiseData <- reactive({
    req(bibliodata())
    df <- as.data.frame(table(bibliodata()$PY))
    colnames(df) <- c("Year", "Freq")
    df$Year <- as.numeric(as.character(df$Year))
    df <- df[order(df$Year), ]
    df$FreqFractnd <- df$Freq / 1000
    df
  })
  
  # UI tab rendering
  output$mainUI <- renderUI({
    if (is.null(bibliodata())) {
      fluidRow(
        box(width = 12, title = "Instructions", status = "info",
            "Please select the input data directory or file and click submit for analysis.")
      )
    } else {
      tabBox(width = 12,
             tabPanel("Summary",
                      box(width = 12, title = "Selected Bibliometric Summary", status = "info",
                          DTOutput("summaryTable"))
             ),
             tabPanel("Data",
                      box(width = 12,
                          fluidRow(
                            column(4,
                                   downloadButton("downloadData", "Download CSV", class = "btn-primary"),
                                   conditionalPanel(
                                     condition = "input.inputMode == 'dir'",
                                     downloadButton("downloadMergedTxt", "Download Merged PubMed File", class = "btn-info")
                                   )
                            )
                          ),
                          div(style = 'overflow-x: auto;', DTOutput("dataTable"))
                      )
             ),
             tabPanel("Year-wise Publication",
                      box(width = 12,
                          fluidRow(
                            column(3,
                                   selectInput("plotType", "Plot Type:", choices = c("Bar", "Line")),
                                   colourpicker::colourInput("plotColor", "Select Color:", value = "#8137b2"),
                                   checkboxInput("showGrid", "Show Grid Lines", value = TRUE),
                                   checkboxInput("swapAxes", "Axis Orientation", value = FALSE),
                                   numericInput("plotWidth", "Plot Width (in):", value = 7, min = 3),
                                   numericInput("plotHeight", "Plot Height (in):", value = 5, min = 3),
                                   numericInput("plotRes", "Resolution (dpi):", value = 300, min = 72),
                                   selectInput("plotFormat", "Download Format:", choices = c("png", "pdf", "tiff")),
                                   downloadButton("downloadPlot", "Download Plot")
                            ),
                            column(9, plotOutput("yearwisePlot", height = "500px"))
                          )
                      )
             ),
             tabPanel("Country-wise Production",
                      box(width = 12, title = NULL, solidHeader = TRUE, status = "primary",
                          collapsible = FALSE,
                          
                          # Collapsible Box 1: Country Production Bar Plot
                          box(title = "Country Production (Top N Countries)",
                              status = "info", width = 12,
                              solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                              fluidRow(
                                column(3,
                                       numericInput("topCountries", "Top N Countries:", value = 10, min = 1),
                                       checkboxInput("showGridCountry", "Show Grid Lines", value = TRUE),
                                       colourpicker::colourInput("scpColor", "SCP Color:", value = "#1b9e77"),
                                       colourpicker::colourInput("mcpColor", "MCP Color:", value = "#d95f02"),
                                       numericInput("countryPlotWidth", "Plot Width (in):", value = 7, min = 3),
                                       numericInput("countryPlotHeight", "Plot Height (in):", value = 5, min = 3),
                                       numericInput("countryPlotRes", "Resolution (dpi):", value = 300, min = 72),
                                       selectInput("countryPlotFormat", "Download Format:", choices = c("png", "pdf", "tiff")),
                                       downloadButton("downloadCountryPlot", "Download Plot")
                                ),
                                column(9,
                                       plotOutput("countryPlot", height = "500px")
                                )
                              )
                          ),
                          
                          # Collapsible Box 2: Country Collaboration Network (to be developed)
                          box(title = "Country Collaboration Network",
                              status = "warning", width = 12,
                              solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                              fluidRow(
                                column(3,
                                       colourpicker::colourInput("nodeColor", "Node Color", value = "#1f77b4"),
                                       colourpicker::colourInput("labelColor", "Label Color", value = "#000000"),
                                       numericInput("labelSize", "Label Size", value = 20, min = 5),
                                       numericInput("edgeWidth", "Edge Width", value = 2, min = 1),
                                       colourpicker::colourInput("edgeColor", "Edge Color", value = "#888888"),
                                       selectInput("layoutType", "Layout", choices = c("layout_nicely", "layout_with_fr", "layout_with_kk", "layout_with_drl")),
                                       downloadButton("downloadNetworkPlot", "Download Network"),
                                       downloadButton("downloadNodeCSV", "Download Node Table (CSV)", class = "btn-success"),
                                       downloadButton("downloadEdgeCSV", "Download Edge Table (CSV)", class = "btn-warning"),
                                       
                                ),
                                column(9,
                                       visNetwork::visNetworkOutput("countryNetworkPlot", height = "600px") %>%
                                         shinycssloaders::withSpinner(color = "#00a65a", type = 6)
                                       
                                )
                              )
                          )
                          
                      )
             )
             
      )
    }
  })
  
  # Summary table
  output$summaryTable <- renderDT({
    req(rv$biblioSummary)
    selected_rows <- c(2:7, 69, 70, 72:74, 76:79)
    datatable(rv$biblioSummary$MainInformationDF[selected_rows, ],
              options = list(dom = 't', ordering = FALSE),
              rownames = FALSE, class = "compact stripe")
  })
  
  # Data table
  output$dataTable <- renderDT({
    req(bibliodata())
    datatable(bibliodata(),
              options = list(pageLength = 10, scrollX = TRUE, scrollY = "400px", paging = TRUE),
              class = 'compact stripe hover',
              rownames = FALSE)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() paste0("bibliometric_data_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(bibliodata(), file, row.names = FALSE)
    }
  )
  
  # Plot output
  output$yearwisePlot <- renderPlot({
    req(yearwiseData())
    df <- yearwiseData()
    df$Year <- as.factor(df$Year)
    
    p <- ggplot(df, aes(x = Year, y = FreqFractnd))
    if (input$plotType == "Bar") {
      p <- p + geom_col(fill = input$plotColor, width = 0.7)
    } else {
      p <- p + geom_line(group = 1, color = input$plotColor, size = 1) +
        geom_point(color = input$plotColor, size = 2)
    }
    
    p <- p + geom_text(aes(label = Freq), vjust = -0.8, size = 3) +
      xlab("Years") + ylab("Frequency (Fractionated)")
    
    base_theme <- theme_minimal(base_size = 13)
    if (!input$showGrid) {
      base_theme <- base_theme +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }
    
    if (input$swapAxes) {
      p <- p + coord_flip()
    }
    
    p + base_theme +
      theme(
        axis.text.x = element_text(angle = if (!input$swapAxes) 90 else 0, hjust = 1),
        axis.text.y = element_text(angle = 0)
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      expand_limits(y = 0.1)
  })
  
  # Download Plot
  output$downloadPlot <- downloadHandler(
    filename = function() paste0("yearwise_plot.", input$plotFormat),
    content = function(file) {
      df <- yearwiseData()
      df$Year <- as.factor(df$Year)
      p <- ggplot(df, aes(x = Year, y = FreqFractnd))
      if (input$plotType == "Bar") {
        p <- p + geom_col(fill = input$plotColor, width = 0.7)
      } else {
        p <- p + geom_line(group = 1, color = input$plotColor, size = 1) +
          geom_point(color = input$plotColor, size = 2)
      }
      p <- p + geom_text(aes(label = Freq), vjust = -0.8, size = 3) +
        xlab("Years") + ylab("Frequency (Fractionated)")
      if (!input$showGrid) {
        p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      }
      if (input$swapAxes) {
        p <- p + coord_flip()
      }
      p <- p + theme_minimal(base_size = 13)
      ggsave(file, plot = p, width = input$plotWidth, height = input$plotHeight,
             dpi = input$plotRes, device = input$plotFormat)
    }
  )
  
  countryData <- reactive({
    req(rv$biblioResults)
    df <- rv$biblioResults$CountryCollaboration
    df <- df[1:min(input$topCountries, nrow(df)), ]
    dfMltd <- reshape2::melt(df[, c("Country", "SCP", "MCP")], id.vars = "Country")
    dfMltd$variable <- factor(dfMltd$variable, levels = rev(sort(unique(dfMltd$variable))))
    dfMltd
  })
  
  output$countryPlot <- renderPlot({
    df <- countryData()
    
    p <- ggplot(df, aes(fill = variable, x = value, y = reorder(Country, value))) +
      geom_bar(position = "stack", stat = "identity") +
      xlab("No of Articles") + ylab("Country") +
      scale_fill_manual(values = c("SCP" = input$scpColor, "MCP" = input$mcpColor)) +
      labs(fill = "Type") +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 0.7, keyheight = 0.7))
    
    base_theme <- theme_minimal(base_size = 12)
    if (!input$showGridCountry) {
      base_theme <- base_theme +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    }
    
    p + base_theme +
      theme(
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
      )
  })
  
  output$downloadCountryPlot <- downloadHandler(
    filename = function() {
      paste0("country_production_plot.", input$countryPlotFormat)
    },
    content = function(file) {
      df <- countryData()
      
      p <- ggplot(df, aes(fill = variable, x = value, y = reorder(Country, value))) +
        geom_bar(position = "stack", stat = "identity") +
        xlab("No of Articles") + ylab("Country") +
        scale_fill_manual(values = c("SCP" = input$scpColor, "MCP" = input$mcpColor)) +
        labs(fill = "Type") +
        guides(fill = guide_legend(reverse = TRUE, keywidth = 0.7, keyheight = 0.7))
      
      base_theme <- theme_minimal(base_size = 12)
      if (!input$showGridCountry) {
        base_theme <- base_theme +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
      }
      
      p <- p + base_theme +
        theme(
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9)
        )
      
      ggsave(file, plot = p, width = input$countryPlotWidth,
             height = input$countryPlotHeight,
             dpi = input$countryPlotRes,
             device = input$countryPlotFormat)
    }
  )
  
  collabGraphData <- reactive({
    # req(rv$biblioResults)
    req(bibliodata())
    
    # data <- rv$biblioResults$Data
    data <- bibliodata()
    if (is.null(data) || !"SR" %in% names(data)) {
      validate("Bibliometric data not fully initialized.")
    }
    
    M <- metaTagExtraction(data, Field = "AU_CO", sep = ";")
    
    if (is.null(M) || !"AU_CO" %in% names(M)) {
      validate("Field 'AU_CO' not found or could not be extracted. Cannot build country collaboration network.")
    }
    
    M$SR <- make.unique(M$SR, sep = "_")
    mat <- biblioNetwork(M, analysis = "collaboration", network = "countries", sep = ";")
    
    g <- igraph::graph_from_adjacency_matrix(as.matrix(mat), mode = "undirected", weighted = TRUE)
    V(g)$degree <- igraph::degree(g)
    g
    print(g)
  })
  
  output$countryNetworkPlot <- renderVisNetwork({
    validate(need(!is.null(collabGraphData()), "Please process bibliometric data first."))
    req(collabGraphData())
    g <- collabGraphData()
    
    nodes <- data.frame(id = igraph::V(g)$name,
                        label = igraph::V(g)$name,
                        value = igraph::V(g)$degree,
                        color = input$nodeColor,
                        font = list(color = input$labelColor, size = input$labelSize),
                        stringsAsFactors = FALSE)
    
    edges <- igraph::as_data_frame(g, what = "edges")
    colnames(edges) <- c("from", "to", "value")
    edges$width <- input$edgeWidth
    edges$color <- input$edgeColor
    
    visNetwork::visNetwork(nodes, edges, height = "600px") %>%
      visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visNetwork::visIgraphLayout(layout = input$layoutType)
  })
  
  # output$downloadNetworkPlot <- downloadHandler(
  #   filename = function() {
  #     paste0("country_collab_network_", Sys.Date(), ".html")
  #   },
  #   content = function(file) {
  #     g <- collabGraphData()
  #     
  #     nodes <- data.frame(id = igraph::V(g)$name,
  #                         label = igraph::V(g)$name,
  #                         value = igraph::V(g)$degree,
  #                         color = input$nodeColor,
  #                         font = list(color = input$labelColor, size = input$labelSize),
  #                         stringsAsFactors = FALSE)
  #     
  #     edges <- igraph::as_data_frame(g, what = "edges")
  #     colnames(edges) <- c("from", "to", "value")
  #     edges$width <- input$edgeWidth
  #     edges$color <- input$edgeColor
  #     
  #     vis <- visNetwork::visNetwork(nodes, edges, height = "600px") %>%
  #       visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  #       visNetwork::visIgraphLayout(layout = input$layoutType)
  #     
  #     visNetwork::visSave(vis, file = file)
  #   }
  # )
  
  output$downloadNetworkPlot <- downloadHandler(
    filename = function() {
      paste0("country_collab_network_", Sys.Date(), ".png")
    },
    content = function(file) {
      g <- collabGraphData()
      
      nodes <- data.frame(
        id = igraph::V(g)$name,
        label = igraph::V(g)$name,
        value = igraph::V(g)$degree,
        color = input$nodeColor,
        font = list(color = input$labelColor, size = input$labelSize),
        stringsAsFactors = FALSE
      )
      
      edges <- igraph::as_data_frame(g, what = "edges")
      colnames(edges) <- c("from", "to", "value")
      edges$width <- input$edgeWidth
      edges$color <- input$edgeColor
      
      vis <- visNetwork::visNetwork(nodes, edges, height = "600px") %>%
        visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
        visNetwork::visIgraphLayout(layout = input$layoutType)
      
      # Save to temporary HTML
      tmp_html <- tempfile(fileext = ".html")
      visNetwork::visSave(vis, file = tmp_html, selfcontained = TRUE)
      
      # Use webshot2 to take snapshot
      webshot2::webshot(url = tmp_html, file = file, selector = "body", vwidth = 1200, vheight = 800)
    }
  )
  
  output$downloadNodeCSV <- downloadHandler(
    filename = function() {
      paste0("network_nodes_", Sys.Date(), ".csv")
    },
    content = function(file) {
      g <- collabGraphData()
      req(g)
      node_df <- data.frame(
        id = igraph::V(g)$name,
        degree = igraph::degree(g)
      )
      write.csv(node_df, file, row.names = FALSE)
    }
  )
  
  output$downloadEdgeCSV <- downloadHandler(
    filename = function() {
      paste0("network_edges_", Sys.Date(), ".csv")
    },
    content = function(file) {
      g <- collabGraphData()
      req(g)
      edge_df <- igraph::as_data_frame(g, what = "edges")
      write.csv(edge_df, file, row.names = FALSE)
    }
  )
  
  output$downloadMergedTxt <- downloadHandler(
    filename = function() {
      paste0("merged_pubmed_", Sys.Date(), ".txt")
    },
    content = function(file) {
      req(mergedPubmedTxtFile())
      file.copy(mergedPubmedTxtFile(), file)
    }
  )
  
  
}


shinyApp(ui, server)
