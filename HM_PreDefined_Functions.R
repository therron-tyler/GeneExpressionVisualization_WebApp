# Version 4.3.1.6_charlie
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
library(ggdendro)
library(stats)
library(RColorBrewer)
library(logger)
library(ragg)

source("./Sourced_Functions/Custom_HierClust_Dendrogram_Functions.R")
source("./Sourced_Functions/Custom_KmeansClust_Function.R")

exp_group_avg_gene_expression_calculator <- function(data, exp_grp) {
  sample_names <- exp_grp$Sample
  groups <- exp_grp$Group
  
  # Filter the data dataframe to include only the columns present in sample_names
  filtered_data <- data %>%
    select(any_of(sample_names))
  
  # Transpose the filtered data so that samples are rows and genes are columns
  transposed_data <- as.data.frame(t(filtered_data))
  
  # Add the information to the transposed data
  transposed_data$Exp_Grp <- groups[match(rownames(transposed_data), sample_names)]
  
  # Group by scler_type and calculate the mean for each gene
  averaged_data <- transposed_data %>%
    group_by(Exp_Grp) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
    
  # Transpose back to the original format with genes as rows and scler_type groups as columns
  averaged_data <- as.data.frame(t(averaged_data))
  
  colnames(averaged_data) <- averaged_data[1,]
  averaged_data <- averaged_data[-1,]
  
  ExperimentalGrp_GE_avg <- cbind(Symbol = data$Symbol, Gene_Symbol = data$Gene_Symbol, averaged_data)
  return(ExperimentalGrp_GE_avg)
}

# ----------------# ----------------# ----------------# ----------------# ----------------
# ----------------# ----------------# ----------------# ----------------# ----------------

groupfile_reformat <- function(group_file) {
  group_df_t <- data.frame(t(group_file))
  organized_hm_groups <- data.frame(Sample = row.names(group_df_t), Group = group_df_t$X1, Color = group_df_t$X2)
  return(organized_hm_groups)
}

# ----------------# ----------------# ----------------# ----------------# ----------------
# ----------------# ----------------# ----------------# ----------------# ----------------


remove_gene_symbol_columns <- function(variable_genes_df) {
  # expression_data <- data.frame(distinct(variable_genes_df), check.names = FALSE)
  expression_data <- variable_genes_df %>% distinct()
  # Set the Gene_Symbol as row names
  rownames(expression_data) <- expression_data$Symbol
  
  # Create a mapping from Ensembl ID to Gene Symbol
  ensembl_to_symbol <- variable_genes_df %>%
    select(Symbol, Gene_Symbol) %>%
    distinct()
  
  # Check for duplicates in Gene_Symbol
  duplicated_symbols <- ensembl_to_symbol$Gene_Symbol[duplicated(ensembl_to_symbol$Gene_Symbol)]
  
  # Append Ensembl IDs to duplicate gene symbols to ensure uniqueness
  ensembl_to_symbol <- ensembl_to_symbol %>%
    mutate(Gene_Symbol = ifelse(Gene_Symbol %in% duplicated_symbols,
                                paste0(Gene_Symbol, "_", Symbol),
                                Gene_Symbol))
  
  # Replace row names in HM_expression_data_with_cluster_info
  new_row_names <- ensembl_to_symbol$Gene_Symbol[match(rownames(expression_data), ensembl_to_symbol$Symbol)]
  
  # Update the row names with the unique gene symbols
  rownames(expression_data) <- new_row_names
  
  # Remove the Gene_Symbol column if it's no longer needed
  expression_data$Gene_Symbol <- NULL
  expression_data$Symbol <- NULL
  
  return(expression_data)
}

# ----------------# ----------------# ----------------# ----------------# ----------------
# ----------------# ----------------# ----------------# ----------------# ----------------

# III
generate_cluster_colors <- function(data, num_clusters) {
  # Determine unique clusters
  unique_clusters <- unique(data$Clusters)
  
  # Ensure you have enough colors for the number of clusters
  if (num_clusters > 24) {
    stop("This method supports up to 24 clusters. Consider revising your clustering method or reducing the number of clusters.")
  }
  
  # Combine colors from multiple palettes to create 24 distinct colors
  cluster_colors <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Accent")[1:(24 - 12 - 8)])
  
  # Map clusters to colors
  cluster_color_map <- setNames(cluster_colors[1:num_clusters], unique_clusters)
  
  # Add the color column to the dataframe
  data$ClusterColor <- cluster_color_map[data$Cluster]
  
  # Arrange data by cluster for consistency
  data <- data %>% arrange(-Cluster)
  
  return(data)
}

# ----------------# ----------------# ----------------# ----------------# ----------------
# ----------------# ----------------# ----------------# ----------------# ----------------

generate_annotations <- function(group_data, exp_color_palette_func) {
  # Extract cell types based on the given pattern
  print("group_data")
  print(group_data)
  # Define colors for groups using the specified color palette function
  exp_colors <- exp_color_palette_func(length(unique(group_data)), begin = 0.5, end = 0.9)
  names(exp_colors) <- unique(group_data)
  col_annotations_exp <- exp_colors[group_data]
  
  return(list(exp_groups = unique(group_data), col_annotations_exp = col_annotations_exp))
}

# ----------------# ----------------# ----------------# ----------------# ----------------
# ----------------# ----------------# ----------------# ----------------# ----------------
capitalize_first_letter <- function(x) {
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

# Function to filter gene expression data based on selected genes
filter_genes <- function(GEtable, selected_genes) {
  # Convert Gene_Symbol to character if it's a factor
  if (is.factor(GEtable$Gene_Symbol)) {
    GEtable$Gene_Symbol <- as.character(GEtable$Gene_Symbol)
  }
  
  # Apply the capitalize_first_letter function safely
  GEtable$Gene_Symbol <- sapply(GEtable$Gene_Symbol, function(x) {
    if (is.na(x) || x == "") {
      return(x)
    } else {
      capitalize_first_letter(x)
    }
  })
  
  # Apply the function to selected_genes as well
  selected_genes <- sapply(selected_genes, function(x) {
    if (is.na(x) || x == "") {
      return(x)
    } else {
      capitalize_first_letter(x)
    }
  })
  
  # Filter the rows where Gene_Symbol matches any of the selected genes
  filtered_GEtable <- GEtable[GEtable$Gene_Symbol %in% selected_genes, ]
  
  return(filtered_GEtable)
}

# ----------------# ----------------# ----------------# ----------------# ----------------
# ----------------# ----------------# ----------------# ----------------# ----------------

generate_space_text <- function(num_columns, base_space = 13) {
  # Adjust the base space depending on the number of columns
  space_multiplier <- floor(num_columns / 10)
  total_space <- base_space + space_multiplier
  return(paste(rep(" ", total_space), collapse = ""))
}

remove_space_text <- function(num_columns, base_space = 10) {
  # Reduce the base space when there are more than 30 columns
  space_multiplier <- floor(num_columns / 6)
  total_space <- max(base_space - space_multiplier, 0)  # Ensure that space count does not go below zero
  return(paste(rep(" ", total_space), collapse = ""))
}


# ----------------# ----------------# ----------------# ----------------# ----------------
# ----------------# ----------------# ----------------# ----------------# ----------------
add_dendrogram_to_HM <- function(clean_variable_genes, hierarchically_clustered_HM) {
  log_info("Entered add_dendrogram_to_HM function.")
  
  # Log the dimensions and structure of inputs
  log_debug("clean_variable_genes dimensions: {dim(clean_variable_genes)}")
  log_debug("hierarchically_clustered_HM class: {class(hierarchically_clustered_HM)}")
  
  hclust_rows <- hierarchical_clustering(clean_variable_genes)
  log_debug("Hierarchical clustering completed.")
  
  log_debug("Dendrogram plot generated.")

  num_rows <- nrow(clean_variable_genes)
  log_debug("Number of rows in clean_variable_genes: {num_rows}")
  
  plot_height <- 800 + 40 * num_rows
  log_debug("Calculated plot height: {plot_height}")
  
  pdf(NULL)
  dendrogram_plotly <- generate_dendrogram(hclust_rows, num_rows, plot_height)
  dev.off()

  # OG solution that works
  dendrogram_plotly <- dendrogram_plotly %>%
    layout(
      autosize = FALSE,
      margin = list(l = 20, r = 0, t = -200, b = 0),
      yaxis = list(range = c(0, num_rows),  # Sync y-axis
                   showline = FALSE,  # Hide the axis line
                   showgrid = FALSE,  # Hide the grid lines
                   zeroline = FALSE,
                   tickvals = seq(-1, num_rows),
                   ticks = "",              # Hide the ticks
                   showticklabels = FALSE), # Hide the zero line
      xaxis = list(showline = FALSE,  # Hide the axis line
                   showgrid = FALSE,  # Hide the grid lines
                   zeroline = FALSE,
                   autorange = "reversed",
                   showticklabels = FALSE)  # Hide the zero line
    )
  
  log_debug("Applied layout to dendrogram plotly object.")

  # Remove row names (y-axis labels) from heatmap
  heatmap_plot <- hierarchically_clustered_HM %>%
    layout(
      yaxis = list(showticklabels = TRUE,
                   tickvals = seq(-1, num_rows, 1),
                   ticks = "")
    )
  
  log_debug("Modified heatmap plot layout.")
  
  combined_plot <- plotly::subplot(
    dendrogram_plotly,
    heatmap_plot,
    nrows = 1,
    widths = c(0.1, 0.9),  # widths for alignment
    margin = 0.02
  )
  log_debug("Combined dendrogram and heatmap into subplot.")
  
  # Apply the config once
  combined_plot <- plotly::config(combined_plot, 
                                  displayModeBar = TRUE, 
                                  modeBarButtonsToRemove = list("toImage"))
  log_info("Final combined plot ready.")
  
  return(combined_plot)
}
