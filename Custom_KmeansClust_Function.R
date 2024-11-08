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

# Function to perform k-means clustering
k_means_clustering <- function(expression_data, num_clusters) {
  # Scale the data before clustering
  
  if (nrow(expression_data) == 0 || ncol(expression_data) == 0) {
    stop("Error: The input expression data is empty. Please provide valid gene expression data.")
  }
  
  # Convert all columns to numeric
  expression_data[] <- lapply(expression_data, function(x) as.numeric(trimws(as.character(x))))
  
  normalized_data_t <- apply(expression_data, 1, min_max_normalization, min_value = -1, max_value = 1)
  
  if (nrow(normalized_data_t) == 0 || ncol(normalized_data_t) == 0) {
    stop("Error: The normalization process returned an empty dataset.")
  }
  
  # Convert the matrix back to a data frame with the same column names
  normalized_data <- t(as.data.frame(normalized_data_t))
  
  num_unique_data_points <- nrow(unique(normalized_data))
  
  if (num_unique_data_points < num_clusters) {
    stop(paste("Error: There are more cluster centers (", num_clusters, 
               ") than distinct data points (", num_unique_data_points, 
               "). Please reduce the number of clusters or broaden your filtering criteria."))
  }
  
  # Perform k-means clustering
  set.seed(40) # For reproducibility
  
  kmeans_result <- kmeans(normalized_data, centers = num_clusters, nstart = 25)
  
  Clusters <- data.frame(Cluster = kmeans_result$cluster)
  
  norm_data_and_clusters <- cbind(normalized_data, Clusters)
  
  return(norm_data_and_clusters)
}
