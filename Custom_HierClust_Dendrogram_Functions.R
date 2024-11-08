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


hierarchical_clustering <- function(GE_data_Symbols) {
  
  # Convert all columns to numeric
  GE_data_Symbols[] <- lapply(GE_data_Symbols, function(x) as.numeric(as.character(x)))
  
  # Check for any NAs introduced by conversion
  if (any(is.na(GE_data_Symbols))) {
    warning("There are NAs in GE_data_Symbols after converting to numeric.")
    print("Rows with NA values:")
    print(which(apply(GE_data_Symbols, 1, function(row) any(is.na(row)))))
  }
  
  # Proceed only if data is numeric
  if (!all(sapply(GE_data_Symbols, is.numeric))) {
    stop("Error: GE_data_Symbols contains non-numeric values after conversion.")
  }
  
  # Convert to matrix
  GE_data_Symbols <- as.matrix(GE_data_Symbols)
  
  # Perform scaling
  scaled_data <- t(scale(t(GE_data_Symbols)))
  
  # Remove rows with NA, NaN, or infinite values
  scaled_data <- scaled_data[!apply(scaled_data, 1, function(row) any(is.na(row) | is.nan(row) | is.infinite(row))), ]
  # Validate that there are enough rows for clustering
  validate(
    need(nrow(scaled_data) >= 2, 
         paste("Error: At least 2 objects are required to perform hierarchical clustering. Currently, there are", nrow(scaled_data), "objects after filtering. Please adjust your filtering criteria."))
  )
  
  # Perform hierarchical clustering
  dist_matrix <- dist(scaled_data)
  set.seed(40)
  hc <- hclust(dist_matrix, method = "complete")
  print("HC process successfully complete")
  return(hc)
}

# -------------

generate_dendrogram <- function(hclust_rows, num_rows, plot_height) {
  # Convert hclust object to a dendrogram
  dend <- stats::as.dendrogram(hclust_rows)
  log_debug("Dendrogram generated.")
  
  # Extract dendrogram segments data
  pdf(NULL)
  dend_data <- ggdendro::dendro_data(dend)
  dev.off()
  log_debug("Dendrogram data converted to dendro_data.")
  
  # Adjust x-values to match the heatmap's y-axis indices
  # Flip the dendrogram to align with the heatmap
  dend_data$segments$x <- dend_data$segments$x - 0.5
  dend_data$segments$xend <- dend_data$segments$xend - 0.5
  log_debug("Dendrogram segments adjusted to align with heatmap indices.")
  
  # Create the dendrogram plot using Plotly
  dendrogram_plotly <- plotly::plot_ly(
    type = 'scatter',
    mode = 'lines',
    showlegend = FALSE,
    hoverinfo = 'none',
    width = 75,
    height = plot_height
  ) %>%
    plotly::add_segments(
      x = dend_data$segments$y,
      y = dend_data$segments$x,
      xend = dend_data$segments$yend,
      yend = dend_data$segments$xend,
      line = list(color = 'black')
    ) 
  
  log_debug("Dendrogram Plotly object generated.")
  return(dendrogram_plotly)
}

