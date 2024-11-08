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
library(RColorBrewer)
library(logger)

source("./Sourced_Functions/PreDefined_Functions.R")
source("./Sourced_Functions/HM_PreDefined_Functions.R")
source("./Sourced_Functions/HM_Plotting_Functions.R")
source("./Sourced_Functions/Custom_KmeansClust_Function.R")

process_heatmap_data <- function(GEtable, grps, dataset_name, clustering_method, cluster_number, dendro, selected_genes, ExpGrp_or_Sample) {
  grps <- groupfile_reformat(grps)
  print(grps)
  
  # Initialize the result list for this dataset
  result <- list(error = FALSE, message = NULL, heatmap = NULL)
  
  # Proceed with GEtable as is (already averaged if necessary)
  filtered_GEtable <- filter_genes(GEtable, selected_genes)
  
  # Get genes that were found in the dataset
  genes_found <- unique(filtered_GEtable$Gene_Symbol)  
  missing_genes <- setdiff(capitalize_first_letter(selected_genes), genes_found)
  
  if (length(genes_found) == 0) {
    # No genes found in this dataset
    result$error <- TRUE
    missing_genes_str <- paste(capitalize_first_letter(selected_genes), collapse = ", ")
    result$message <- paste0(
      "None of the selected genes are present in dataset: '", dataset_name, "'.",
      " Missing genes: ", missing_genes_str, "."
    )
    return(result)
  } else if (length(missing_genes) > 0) {
    # Some genes are missing
    missing_genes_str <- paste(missing_genes, collapse = ", ")
    result$message <- paste0(
      "The following genes are not expressed: ",
      missing_genes_str)
    # Proceed to generate the heatmap with the available genes
  }
  
  # Remove symbol columns from gene expression table
  cleaned_data <- remove_gene_symbol_columns(filtered_GEtable)
  
  # Determine clustering method
  if (clustering_method == "kmeans") {
    num_unique_data_points <- nrow(unique(cleaned_data))
    
    if (num_unique_data_points <= cluster_number) {
      result$error <- TRUE
      result$message <- paste0(
        "The number of k-means clusters (", cluster_number, 
        ") must be less than the number of genes (", num_unique_data_points, ") in dataset: '", dataset_name, "'."
      )
      return(result)
    }
    
    # Proceed with K-Means clustering
    kmeans_clusters_and_variableGE_data <- k_means_clustering(cleaned_data, cluster_number)
    HM_expression_data_with_cluster_info <- generate_cluster_colors(kmeans_clusters_and_variableGE_data, cluster_number)
    heatmap <- kmeans_heatmap(HM_expression_data_with_cluster_info, grps, ExpGrp_or_Sample)
    result$heatmap <- heatmap
    
  } else if (clustering_method == "hierarchical") {
    num_data_points <- nrow(cleaned_data)
    
    if (num_data_points < 2) {
      result$error <- TRUE
      result$message <- paste0(
        "At least two genes are required for hierarchical clustering in dataset: '", dataset_name, "'."
      )
      return(result)
    }
    
    heatmap <- HierClust_HM(cleaned_data, grps, ExpGrp_or_Sample)
    
    # Proceed with hierarchical clustering
    if (dendro == "Yes") {
      hm_no_dendro <- heatmap
      heatmap <- add_dendrogram_to_HM(cleaned_data, hm_no_dendro)
    } 
    
    result$heatmap <- heatmap
    
  } else {
    # No clustering
    heatmap <- no_clustering_heatmap(cleaned_data, grps, ExpGrp_or_Sample)
    result$heatmap <- heatmap
  }
  
  return(result)
}

