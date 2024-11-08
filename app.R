# Gene Visualization: Winter Lab Database of bulk RNA-seq Data 11/1/2024
# Tyler Therron
# Version 4.3.1.6_charlie

#  ======================================. Libraries ==================

library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(readr)
library(shiny.fluent)
library(htmlwidgets) 
library(webshot) 
library(future)
library(DT)
library(future.apply)
library(RColorBrewer)
library(logger)
library(ragg)

options(shiny.maxRequestSize = 1000 * 1024^2)
#  ======================================. Global Section for help functions and data ==================

source("./Sourced_Functions/PreDefined_Functions.R")
source("./Sourced_Functions/HM_PreDefined_Functions.R")
source("./Sourced_Functions/HM_Plotting_Functions.R")
source("./Sourced_Functions/Heatmap_Integration_Function.R")

#  ======================================. Define UI ==================
# Configure logger
# log_file_path <- "/Users/ttm3567/Documents/October2024/GeneExpressionVisualization_v4.3.1_beta_logging/internal_server_GEV_app.log"
log_file_path <- "/srv/shiny-server/internal_server_GEV_app.log"
if (!dir.exists(dirname(log_file_path))) {
  dir.create(dirname(log_file_path), recursive = TRUE)
}
log_appender(appender_file(log_file_path))
log_threshold(DEBUG)
log_layout(layout_glue_generator(format = '{time} {level} [{namespace}] {fn}(): {msg}'))

ui <- fluidPage(
  # Include custom styles in the head section
  tags$head(
    tags$style(HTML("
            .grid-container {
              display: grid;
              grid-template-columns: repeat(2, 1fr);
              grid-gap: 20px;
              padding: 20px;
            }
            .grid-item {
              background-color: #f8f9fa;  /* Light grey background */
              padding: 10px;
              border-radius: 5px;
              box-shadow: 0 2px 4px rgba(0,0,0,0.1);  /* Subtle shadow for depth */
              overflow: hidden; /* Ensures the content fits in the grid area */
            }
        "))
  ),
  titlePanel("Winter Lab Database of Bulk RNA-seq"),
  tabsetPanel(
    id = "main_tabs",
    tabPanel(
      "Introduction",
      h2("Welcome to the Winter Lab Database of Bulk RNA-seq!"),
      br(),
      p(
        style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
        "This application contains published and unpublished Bulk RNA-seq datasets from mice.
      The primary purpose of this tool is to explore normalized gene expression levels among experimental groups from preloaded or user-loaded Bulk RNA-seq datasets."
      ),
      br(),
      br(),
      fluidRow(tags$div(
        style = "display: flex; justify-content: center; margin-left: 40px; margin-right: 40px;",
        DT::dataTableOutput("Dataset_Info")
      )),
      br(),
      br(),
      fluidRow(
        p(
          style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
          "For instructions on how to properly format datatables and create groupfiles, click the link below:"
        ),
        tags$a(
          style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
          id = "newapp",
          href = "https://winterlab-webapps.fsm.northwestern.edu/Template_GroupFileMaker_v3_App_20231120/",
          "Template Data and Experimental Group Formatting App."
        ),
        tags$script(
          HTML(
            "$(document).on('click', '#newapp', function (e) {
                      e.preventDefault();
                      window.open($(this).attr('href'));
                  });"
          )
        ),
        br(),
        br(),
        br(),
        p(
          style = "font-size: 20px; margin-left: 40px; margin-right: 40px;",
          "Below is a table with definitions of all the acronyms used in this web application."
        ),
        tags$div(
          style = "display: flex; justify-content: left; margin-left: 40px; margin-right: 40px;",
          tableOutput("legend_dataset")
        ),
        br(),
        br()
      )
    ),
    tabPanel(
      "Data Visualization",
      br(),
      sidebarLayout(
        sidebarPanel(
          div(
            p(
              "To visualize gene expression in mouse datasets.\n Select one or two gene(s) of interest, and \nplot an interactive boxplot for gene expression data"
            )
          ),
          checkboxInput("select_all", "Select All Preloaded Data"),
          checkboxGroupInput("gene_expression_data", "Preloaded Data:", choices = NULL),
          checkboxInput("uploaded_data", "Uploaded Data:"),
          div(
            p(
              "To download a plot as a PNG image, click the camera icon in the menu at the top-right of the rendered plot.",
              style = "font-style: italic; font-size: smaller;"
            )
          ),
          conditionalPanel(
            condition = "input.uploaded_data == true",
            fileInput(
              'uploaded_data_file1', 'Choose CSV File',
              accept = c(
                'text/csv',
                'text/comma-separated-values,text/plain',
                '.csv'
              )
            ),
            # textInput("uploaded_data_title", "Title for First Dataset:"),
            fileInput(
              "groups_File1", "Upload sample groups",
              multiple = FALSE,
              accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
            )
          ),
          conditionalPanel(
            condition = "input.uploaded_data == true",
            tags$a(id = "newapp", href = "https://winterlab-webapps.fsm.northwestern.edu/Template_GroupFileMaker_v3_App_20231120/", "Create Formatted Experimental Group File."),
            tags$script(
              HTML(
                "$(document).on('click', '#newapp', function (e) {
                      e.preventDefault();
                      window.open($(this).attr('href'));
                  });"
              )
            )
          ),
          radioButtons("cpm_or_fpkm", "CPM or FPKM preloaded data repository:", c("CPM" = "cpm", "FPKM" = "fpkm"), selected = "cpm"),
          radioButtons("colSel_pri", "Select Gene Nomenclature:", c("Ensembl ID" = "one", "Gene Symbol" = "two"), selected = "two"),
          radioButtons("gene_select_method", "Upload a gene list or select individual genes to visualize:", c("Individual Genes" = "gene", "Gene List CSV" = "csv"), selected = "gene"),
          conditionalPanel(
            condition = "input.gene_select_method == 'csv'",
          fileInput("gene_file", "Upload CSV of Genes", 
                    accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))),
          conditionalPanel(
            condition = "input.gene_select_method == 'gene'",
          uiOutput('select_pri')),
          width = 3
        ),
        mainPanel(
          uiOutput("dataset_tabs")
        )
      )
    ),
    tabPanel( "Gene Expression Heatmaps",
              sidebarLayout(
                sidebarPanel(
                  selectInput("cluster_hm_method", "Select Method to Cluster Heatmap Genes by:", choices = c("K-Means" = "kmeans",
                                                                                                             "Hierarchical" = "hierarchical",
                                                                                                             "No Clustering" = "none"), selected = "kmeans"),
                  conditionalPanel( # CPM files - this is where the groupmaker app link will be
                    condition = "input.cluster_hm_method == 'kmeans'",
                    numericInput("number_o_clusters", "Number of K-Means Clusters:", value = 2, min = 0, step = 1)),
                  conditionalPanel( # CPM files - this is where the groupmaker app link will be
                    condition = "input.cluster_hm_method == 'hierarchical'",
                    selectInput("dendrogram", "Dendrogram?", choices = c("No", "Yes"), selected = "No")),
                  selectInput("ExpGrp_or_Sample", "Visualize by Experimental Groups or Samples?", choices = c("Group", "Sample"), selected = "Group"),
                  width = 3),
                mainPanel(br(),
                          tags$div(
                            style = "overflow-y: scroll; height: 1200px;",
                          uiOutput("heatmap_tabs")), br(),
                          div(
                            p(
                              "(The heatmap is in a scrollable window.)",
                              style = "font-style: italic; font-size: smaller;"
                            )
                          )
                          , width = 9)))
  )
) # UI script
# ======= SERVER ====================================================

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  log_info("Session started: {session$token}")
  log_info("R version: {R.Version()$version.string}")
  log_info("ggdendro version: {packageVersion('ggdendro')}")
  log_info("ggplot2 version: {packageVersion('ggplot2')}")
  log_info("dendextend version: {packageVersion('dendextend')}")
  
  selected_genes <- reactiveValues(genes = NULL)
  
  gene_file_genes <- reactive({
    req(input$gene_file)
    
    # Read the uploaded file and extract the gene symbols
    gene_data <- read.csv(input$gene_file$datapath, header = TRUE, stringsAsFactors = FALSE)
    
    # Assuming the gene symbols are in the first column
    gene_vector <- tolower(gene_data[, 1])  # Ensure it is a character vector and lowercase
    
    # Return the vector of gene symbols
    return(gene_vector)
  })
  
  processed_datasets <- reactiveValues(data = list(), groups = list())
  
  # parallelized data read 
  observeEvent(input$gene_expression_data, {
    log_debug("The dataset checkboxes have been clicked by user. {session$token}: {input$gene_expression_data}")
    selected_data <- input$gene_expression_data
    new_data <- setdiff(selected_data, names(processed_datasets$data))
    
    if (length(new_data) > 0) {
      withProgress(message = 'Loading datasets...', value = 0, {
        future_lapply(seq_along(new_data), function(i) {
          dataset <- new_data[i]
          path <- if (input$cpm_or_fpkm == "cpm") {
            paste0('./Primary/', dataset, '/Data.csv')
          } else {
            paste0('./FPKM_expression_files/', dataset, '/Data.csv')
          }
          
          validate(need(file.exists(path), "Data file does not exist."))
          
          data <- read.csv(path, check.names = FALSE)
          out_filtered_genes <- data[apply(Filter(is.numeric, data), 1, function(row) any(row != 0)), ]
          processed_datasets$data[[dataset]] <- out_filtered_genes
          
          group_path <- if (input$cpm_or_fpkm == "cpm") {
            paste0('./Primary/', dataset, '/GroupFile.csv')
          } else {
            paste0('./FPKM_expression_files/', dataset, '/GroupFile.csv')
          }
          
          validate(need(file.exists(group_path), "Group file does not exist."))
          
          group_data <- read.csv(group_path)
          processed_datasets$groups[[dataset]] <- group_data
          
          # Update progress bar
          incProgress(1/length(new_data), detail = paste("Loaded", dataset))
        }, future.seed = TRUE)
      })
    }
  })
  
  observe({
    if (input$uploaded_data == TRUE) {
      showNotification("Uploaded data checkbox checked", type = "message")
      
      # Ensure required inputs are provided
      req(input$uploaded_data_file1, input$groups_File1)
      
      data_path <- input$uploaded_data_file1$datapath
      group_path <- input$groups_File1$datapath
      dataset_name <- "User Uploaded Data"
      
      # Check if files exist at the provided paths
      validate(
        need(file.exists(data_path), "User data file does not exist."),
        need(file.exists(group_path), "User group file does not exist.")
      )
      
      showNotification("Files found and being read", type = "message")
      
      # Read the data file
      data <- tryCatch(
        {
          read.csv(data_path, check.names = FALSE)
        },
        error = function(e) {
          showNotification(paste("Error reading data file:", e$message), type = "error")
          print(paste("Error reading data file:", e$message))
          return(NULL)
        }
      )
      
      # Check if data was read successfully
      if (!is.null(data)) {
        
        showNotification(paste("Data:", nrow(data), "rows read"), type = "message")
        
        # parallelized Filter and store the data
        out_filtered_genes <- data[future_apply(Filter(is.numeric, data), 1, function(row) any(row != 0)), ]
        processed_datasets$data[[dataset_name]] <- out_filtered_genes
        
        showNotification("Data processed and filtered", type = "message")
        
        # Read the group file
        group_data <- tryCatch(
          {
            read.csv(group_path)
          },
          error = function(e) {
            showNotification(paste("Error reading group file:", e$message), type = "error")
            print(paste("Error reading group file:", e$message))
            return(NULL)
          }
        )
        
        # Check if group data was read successfully
        if (!is.null(group_data)) {
          
          showNotification(paste("Group data:", nrow(group_data), "rows read"), type = "message")
          
          # Store the group data
          processed_datasets$groups[[dataset_name]] <- group_data
          
          showNotification("Group data processed and stored", type = "message")
          
        } else {
          print("Group data is NULL")
        }
      } else {
        print("Data is NULL")
      }
    } else {
      print("upload data box not checked")
    }
  })
  
  # change to include user uploaded data as well
  preloaded_data_readin_fxn <- reactive({
    
    all_data <- processed_datasets$data
    all_groups <- processed_datasets$groups
    
    return(list(data = all_data, groups = all_groups))
  })
  
  
  datasetSelected <- reactiveVal(NULL)
  
  lastDataset <- reactiveValues(data = NULL)
  
  observe({
    lastDataset$data <- input$gene_expression_data
  })
  
  observeEvent(input$cpm_or_fpkm, {
    datasets <- datasets_info("./FPKM_expression_files", "./Primary")
    
    if (input$cpm_or_fpkm == "cpm") {
      cpm_path <- "./Primary"
      datasetSelected("cpm")
      cpm_choice_list <- c(sort_numeric(list.files(path = cpm_path, full.names = FALSE, recursive = FALSE)))
      updateCheckboxGroupInput(session, "gene_expression_data", choices = cpm_choice_list, selected = if (is.null(lastDataset$data) || length(lastDataset$data) == 0) { cpm_choice_list[1] } else { lastDataset$data })
    } else if (input$cpm_or_fpkm == "fpkm") {
      fpkm_path <- "./FPKM_expression_files"
      datasetSelected("fpkm")
      fpkm_choice_list <- c(sort_numeric(list.files(path = fpkm_path, full.names = FALSE, recursive = FALSE)))
      updateCheckboxGroupInput(session, "gene_expression_data", choices = fpkm_choice_list, selected = if (is.null(lastDataset$data) || length(lastDataset$data) == 0) { fpkm_choice_list[1] } else { lastDataset$data })
    }
  })
  
  # Handle the "Select All" functionality
  observeEvent(input$select_all, {
    req(input$cpm_or_fpkm)
    cpm_path <- "./Primary"
    fpkm_path <- "./FPKM_expression_files"
    cpm_choice_list <- c(sort_numeric(list.files(path = cpm_path, full.names = FALSE, recursive = FALSE)))
    fpkm_choice_list <- c(sort_numeric(list.files(path = fpkm_path, full.names = FALSE, recursive = FALSE)))
    if (input$cpm_or_fpkm == "cpm") {
      choices <- cpm_choice_list
    } else if (input$cpm_or_fpkm == "fpkm") {
      choices <- fpkm_choice_list
    } else {
      choices <- character(0)
    }
    
    if (input$select_all == TRUE) {
      updateCheckboxGroupInput(session, "gene_expression_data", selected = choices)
    } else {
      updateCheckboxGroupInput(session, "gene_expression_data", selected = character(0))
    }
  })
  
  observeEvent(input$gene_expression_data, {
    req(input$cpm_or_fpkm)
    cpm_path <- "./Primary"
    fpkm_path <- "./FPKM_expression_files"
    cpm_choice_list <- c(sort_numeric(list.files(path = cpm_path, full.names = FALSE, recursive = FALSE)))
    fpkm_choice_list <- c(sort_numeric(list.files(path = fpkm_path, full.names = FALSE, recursive = FALSE)))
    
    selected_datasets <- input$gene_expression_data
    if (input$cpm_or_fpkm == "cpm") {
      choices <- cpm_choice_list
    } else if (input$cpm_or_fpkm == "fpkm") {
      choices <- fpkm_choice_list
    } else {
      choices <- character(0)
    }
    
    if (length(selected_datasets) == length(choices)) {
      updateCheckboxInput(session, "select_all", value = TRUE)
    } else {
      updateCheckboxInput(session, "select_all", value = FALSE)
    }
  })
  
  observe({
    datasets <- datasets_info("./FPKM_expression_files", "./Primary")
    if (!is.null(lastDataset$data) && length(lastDataset$data) > 0 && any(lastDataset$data %in% datasets$unavailable)) {
      updateRadioButtons(session, "cpm_or_fpkm", selected = "cpm", choices = c(CPM = "cpm"))
    } else {
      updateRadioButtons(session, "cpm_or_fpkm", selected = input$cpm_or_fpkm, choices = c(CPM = "cpm", FPKM = "fpkm"))
    }
  })
  
  observeEvent(input$cpm_or_fpkm, {
    if (input$cpm_or_fpkm == "cpm") {
      datasetSelected("cpm")
    } else if (input$cpm_or_fpkm == "fpkm") {
      datasetSelected("fpkm")
    }
  })
  
  selNum_pri <- reactive(switch(input$colSel_pri, one = 1, two = 2))
  
  title_pri <- reactive({
    get_data_title(input$gene_expression_data, input$uploaded_data_title)
  })
  
  all_ids_pri <- reactive({
    preloaded_data <- preloaded_data_readin_fxn()
    all_data <- preloaded_data$data
    
    unique_ids <- unique(unlist(future_lapply(all_data, function(dataset) {
      unique(dataset[, selNum_pri()])
    })))
    
    return(unique_ids)
  })
  
  choiceList_pri <- reactive(tolower(all_ids_pri()))
  
  preprocessing_selected_genes_datasets <- reactive({
    combo_data <- preloaded_data_readin_fxn()

    data_in <- combo_data$data

    selected_genes <- selected_genes$genes

    results_list <- future_lapply(data_in, function(dataset) {
      available_genes <- tolower(dataset[, selNum_pri()])
      validate(need(any(available_genes %in% selected_genes), "None of the currently selected genes are available in the data."))

      future_lapply(selected_genes, function(gene) {
        gene_data <- dataset[tolower(dataset[[selNum_pri()]]) == tolower(gene), ]
        process_gene_data(gene_data, gene, selNum_pri)
      })
    })

    # Use assign_genes_to_sublist to name the gene lists
    named_results_list <- assign_names_to_sublist(results_list, selected_genes)
    return(named_results_list)
  })
  
  
  formattedData_pri <- reactive({
    results_list <- preprocessing_selected_genes_datasets()
    combo_data <- preloaded_data_readin_fxn()
    groups <- combo_data$groups

    group_and_data <- future_mapply(function(dataset_genes_list, groupfile_list) {
      format_and_merge_UI_chosen_genes(dataset_genes_list, groupfile_list)
    }, results_list, groups, SIMPLIFY = FALSE)

    return(group_and_data)
  })
  
  
  cpm_fpkm_labeler_Yaxis <- reactive({
    if (input$cpm_or_fpkm == "cpm") {
      label <- "CPM"
    } else if (input$cpm_or_fpkm == "fpkm") {
      label <- "FPKM"
    }
    return(label)
  })
  
  observe({
    if (input$gene_select_method == "gene") {
      selected_genes$genes <- input$mygene_pri
    } else {
      selected_genes$genes <- gene_file_genes()
    }

  })
  
  
  # GOOD - original gene choice
  output$select_pri = renderUI({
    req(preloaded_data_readin_fxn())
    choiceList_pri <- choiceList_pri()

    
    selectizeInput("mygene_pri",
                   label = h5("Select Gene of Interest to Visualize:"),
                   multiple = TRUE,
                   selected = if (is.null(input$mygene_pri) || length(input$mygene_pri) == 0) { "" }
                              else { selected_genes$genes },
                   choices = choiceList_pri,
                   options = list(create = TRUE, placeholder = 'Type here or choose from the list'))
  })
  
  
  # parallelized
  output$dataset_tabs <- renderUI({
    req(input$gene_expression_data)
    
    # Initialize the list of datasets to include the preloaded datasets
    datasets <- input$gene_expression_data
    
    # Include the uploaded dataset if the checkbox is checked
    if (input$uploaded_data == TRUE) {
      datasets <- c(datasets, "User Uploaded Data")
    }
    
    # Get selected genes
    selected_genes <- selected_genes$genes  # Ensure this is accessible
    
    # Get the datasets
    combo_data <- preloaded_data_readin_fxn()
    GEtable_list <- combo_data$data
    
    # Subset GEtable_list to only include selected datasets
    GEtable_list <- GEtable_list[datasets]
    
    # Generate tabs for each dataset
    dataset_tabs <- lapply(datasets, function(dataset) {
      # Get the data for this dataset
      GEtable <- GEtable_list[[dataset]]
      
      filtered_GEtable <- filter_genes(GEtable, selected_genes)
      
      # Get the genes present in the dataset
      genes_found <- unique(filtered_GEtable$Gene_Symbol)
      
      # Find missing genes (case-insensitive comparison)
      missing_genes <- setdiff(capitalize_first_letter(selected_genes), genes_found)
      
      # Create the message if there are missing genes
      if (length(missing_genes) > 0) {
        missing_genes_message <- paste(
          "The following genes could not be visualized in this dataset: ",
          paste(missing_genes, collapse = ", ")
        )
        missing_genes_ui <- div(
          class = "missing-genes-message",
          p(missing_genes_message)
        )
      } else {
        missing_genes_ui <- NULL
      }
      
      # Create the tab panel for the dataset
      tabPanel(
        title = dataset,
        # Include the missing genes message if any
        missing_genes_ui,
        # The gene plots UI output
        uiOutput(paste0("gene_plots_", gsub("[^a-zA-Z0-9]", "_", dataset)))
      )
    })
    
    # Return the tabset panel containing all dataset tabs
    do.call(tabsetPanel, dataset_tabs)
  })
  
  formatted_data <- reactiveValues(data = NULL)
  
  observe({
    formatted_data$data <- formattedData_pri()
  })
  
  # parallelized
  observe({
    req(input$gene_expression_data, datasetSelected())
    
    # Initialize the list of datasets to include the preloaded datasets
    datasets <- input$gene_expression_data
    
    # Include the uploaded dataset if the checkbox is checked
    if (input$uploaded_data) {
      datasets <- c(datasets, "User Uploaded Data")
    }
    
    new_data <- setdiff(datasets, names(formatted_data$data))
    
    if (length(new_data) > 0) {
      results_list <- preprocessing_selected_genes_datasets()
      combo_data <- preloaded_data_readin_fxn()
      groups <- combo_data$groups
      
      # Use future_mapply to parallelize the processing
      group_and_data <- future_mapply(function(dataset_genes_list, groupfile_list) {
        format_and_merge_UI_chosen_genes(dataset_genes_list, groupfile_list)
      }, results_list, groups, SIMPLIFY = FALSE)
      
      formatted_data <<- c(formatted_data, group_and_data)
    }
  })
  
  observeEvent(formattedData_pri(), {
    formatted_data <- formattedData_pri()
    
    withProgress(message = 'Rendering plots...', value = 0, {
      total_datasets <- length(names(formatted_data))
      dataset_count <- 0
      
      # First loop: Create UI outputs for plots
      future_lapply(names(formatted_data), function(dataset_name) {
        output[[paste0("gene_plots_", gsub("[^a-zA-Z0-9]", "_", dataset_name))]] <- renderUI({
          gene_plots <- lapply(names(formatted_data[[dataset_name]]), function(gene_name) {
            plot_data <- formatted_data[[dataset_name]][[gene_name]]
            
            # Check for the 'Count' column
            if (!("Count" %in% names(plot_data))) {
              # Skip this gene
              return(NULL)
            }
            
            plot_output_id <- paste0(
              "plotly_plot_",
              gsub("[^a-zA-Z0-9]", "_", dataset_name),
              "_",
              gsub("[^a-zA-Z0-9]", "_", gene_name)
            )
            
            div(
              plotlyOutput(outputId = plot_output_id, height = "800px", width = "100%"),
              class = "grid-item"
            )
          })
          
          # Remove NULLs
          gene_plots <- Filter(Negate(is.null), gene_plots)
          
          div(
            class = "grid-container",
            do.call(tagList, gene_plots)
          )
        })
        
        dataset_count <- dataset_count + 1
        incProgress(1 / total_datasets, detail = paste("Preparing plots for", dataset_name))
      }, future.seed = TRUE)
      
      # Second loop: Render plots
      future_lapply(names(formatted_data), function(dataset_name) {
        lapply(names(formatted_data[[dataset_name]]), function(gene_name) {
          plot_data <- formatted_data[[dataset_name]][[gene_name]]
          
          # Check for the 'Count' column
          if (!("Count" %in% names(plot_data))) {
            # Skip this gene
            return(NULL)
          }
          
          plot_output_id <- paste0(
            "plotly_plot_",
            gsub("[^a-zA-Z0-9]", "_", dataset_name),
            "_",
            gsub("[^a-zA-Z0-9]", "_", gene_name)
          )
          
          output[[plot_output_id]] <- renderPlotly({
            plot_the_data(plot_data, dataset_name, cpm_fpkm_labeler_Yaxis(), gene_name)
          })
        })
        
        dataset_count <- dataset_count + 1
        incProgress(1 / total_datasets, detail = paste("Rendering plots for", dataset_name))
      }, future.seed = TRUE)
    })
  })
  
  output$Dataset_Info <- renderDataTable({
    datatable(
      dataset_info,
      rownames = FALSE,
      escape = FALSE,
      options = list(
        autoWidth = FALSE,
        columnDefs = list(list(targets = "_all", className = "dt-center")),
        pageLength = -1,
        dom = 't'
      )
    )
  }, server = FALSE)
  
  output$legend_dataset <- renderTable(legend_info, bordered = TRUE, striped = TRUE, rownames = FALSE)
  
  # HM tab reactive functions ---------------------
  
  ## REACTIVE heatmaps
  averaged_data_list <- reactive({
    # Create a progress bar
    withProgress(message = "Averaging data...", value = 0, {
      # Get the original data and groups
      combo_data <- preloaded_data_readin_fxn()
      GEtable_list <- combo_data$data
      grps_list <- combo_data$groups
      selected_genes <- selected_genes$genes
      
      total_datasets <- length(GEtable_list)
      progress_increment <- 1 / total_datasets  # Calculate the increment for each dataset
      
      # Initialize an empty list to store the results
      averaged_data <- list()
      
      # Loop over the datasets, updating the progress bar
      for (i in seq_along(GEtable_list)) {
        dataset_name <- names(GEtable_list)[i]
        GEtable <- GEtable_list[[i]]
        grps <- grps_list[[i]]
        
        tryCatch({
          # Process the data
          grps <- groupfile_reformat(grps)
          filtered_GEtable <- filter_genes(GEtable, selected_genes)
          GEtable_avg <- exp_group_avg_gene_expression_calculator(filtered_GEtable, grps)
          
          # Store the result
          averaged_data[[dataset_name]] <- GEtable_avg
        }, error = function(e) {
          # Handle the error (e.g., log it or show a message)
          showNotification(paste("Error processing", dataset_name, ":", e$message), type = "error")
        })
        
        # Update the progress bar
        incProgress(progress_increment, detail = paste("Processing", dataset_name))
      }
      
      return(averaged_data)
    })
  })
  

  combined_HM <- reactive({
    if (input$main_tabs != "Gene Expression Heatmaps") {
      return(NULL)  # Exit the reactive function if not on Heatmap tab
    }
    
    clustering_method <- input$cluster_hm_method
    cluster_number <- input$number_o_clusters
    dendro <- input$dendrogram
    combo_data <- preloaded_data_readin_fxn()
    GEtable_list <- combo_data$data
    grps_list <- combo_data$groups
    ExpGrp_or_Sample <- input$ExpGrp_or_Sample
    selected_genes <- selected_genes$genes
    
    # Use averaged data if "Group" is selected
    if (ExpGrp_or_Sample == "Group") {
      data_list <- averaged_data_list()
      if (is.null(data_list)) {
        # Data not yet averaged or "Group" not selected
        return(NULL)
      }
    } else {
      data_list <- GEtable_list
    }
    
    # Use Map to apply the helper function over data_list and grps_list
    heatmap_list <- Map(process_heatmap_data, 
                        data_list, grps_list, names(data_list), 
                        MoreArgs = list(clustering_method = clustering_method, 
                                        cluster_number = cluster_number, 
                                        dendro = dendro, 
                                        selected_genes = selected_genes, 
                                        ExpGrp_or_Sample = ExpGrp_or_Sample))
    
    return(heatmap_list)
  })
  
  
  # Generate the heatmap tabs dynamically using combined_HM()
  output$heatmap_tabs <- renderUI({
    req(input$gene_expression_data)
    
    # Initialize the list of datasets to include the preloaded datasets
    datasets <- input$gene_expression_data
    
    # Include the uploaded dataset if the checkbox is checked
    if (input$uploaded_data == TRUE) {
      datasets <- c(datasets, "User Uploaded Data")
    }
    
    # Create tabs for each dataset
    heatmap_tabs <- lapply(datasets, function(dataset) {
      tabPanel(
        title = dataset,
        plotlyOutput(outputId = paste0("heatmap_plot_", gsub("[^a-zA-Z0-9]", "_", dataset)), height = "600px")
      )
    })
    
    # Combine all the tabs into one tabset panel
    do.call(tabsetPanel, heatmap_tabs)
  })
  
  # Render heatmap plots using combined_HM()
  observeEvent(combined_HM(), {
    heatmap_list <- combined_HM()
    
    lapply(seq_along(heatmap_list), function(i) {
      dataset_name <- names(heatmap_list)[i]
      result <- heatmap_list[[i]]
      
      plot_output_id <- paste0("heatmap_plot_", gsub("[^a-zA-Z0-9]", "_", dataset_name))
      log_debug("Heatmap generated. {session$token}: {plot_output_id}")
      
      # Render each heatmap in its respective tab
      output[[plot_output_id]] <- renderPlotly({
        req(result)  # Ensure result is available
        
        if (result$error) {
          # Create an empty plot and display the error message as an annotation
          plot_ly() %>%
            layout(
              annotations = list(
                list(
                  x = 0.5,
                  y = 0.5,
                  text = result$message,
                  showarrow = FALSE,
                  xref = "paper",
                  yref = "paper",
                  xanchor = "center",
                  yanchor = "middle",
                  align = "center",
                  font = list(size = 16, color = "grey20", family = "Arial Italic")
                )
              ),
              xaxis = list(showticklabels = FALSE, zeroline = FALSE, showline = FALSE, showgrid = FALSE),
              yaxis = list(showticklabels = FALSE, zeroline = FALSE, showline = FALSE, showgrid = FALSE)
            )
        } else {
          # If there's a message (e.g., some genes are missing), include it as an annotation below the heatmap
          if (!is.null(result$message)) {
            heatmap <- result$heatmap %>%
              layout(
                annotations = list(
                  list(
                    x = 0.5,
                    y = -0.25,  # Position below the plot
                    text = result$message,
                    showarrow = FALSE,
                    xref = "paper",
                    yref = "paper",
                    xanchor = "center",
                    yanchor = "top",
                    align = "center",
                    font = list(size = 14, color = "grey40", family = "Arial Italic"),
                    valign = "top"
                  )
                ),
                margin = list(b = 200)  # Increase bottom margin to accommodate the annotation
              )
            heatmap
          } else {
            result$heatmap
          }
        }
      })
    })
  })
  
  session$onSessionEnded(function() {
    log_info("Session ended: {session$token}")
  })
} # SERVER script bracket

# Run the application 
shinyApp(ui = ui, server = server)
