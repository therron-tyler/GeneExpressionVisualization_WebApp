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

source("./Sourced_Functions/HM_PreDefined_Functions.R")

prepare_heatmap_data <- function(data, grps, Group_or_Sample, exclude_cols = NULL) {
  # Remove specified columns if needed
  if (!is.null(exclude_cols)) {
    data <- data %>% select(-one_of(exclude_cols))
  }
  
  # Determine col_ID based on Group_or_Sample
  if (Group_or_Sample == "Group") {
    col_ID <- unique(grps$Group)
  } else {
    col_ID <- grps$Sample
  }
  
  # Order columns
  column_order <- match(col_ID, colnames(data))
  data_ordered <- data[, column_order, drop = FALSE]
  
  # Convert columns to numeric
  data_ordered[] <- lapply(data_ordered, function(x) as.numeric(trimws(as.character(x))))
  
  return(data_ordered)
}

generate_heatmap_annotations <- function(grps, Group_or_Sample) {
  if (Group_or_Sample == "Group") {
    column_annotation_input <- unique(grps$Group)
  } else {
    column_annotation_input <- grps$Group
  }
  
  annotations <- generate_annotations(column_annotation_input, viridis::turbo)
  EXP_groups <- annotations$exp_groups
  col_annotations <- annotations$col_annotations_exp
  
  return(list(EXP_groups = EXP_groups, col_annotations = col_annotations))
}

calculate_heatmap_parameters <- function(num_genes, num_columns) {
  if (num_columns <= 30) {
    dynamic_space <- generate_space_text(num_columns)
  } else {
    dynamic_space <- remove_space_text(num_columns)
  }
  
  plot_height <- 800 + 40 * num_genes 
  
  font_size <- ifelse(num_genes > 100, 6, 
                      ifelse(num_genes > 50, 8, 
                             ifelse(num_genes > 25, 10, 
                                    ifelse(num_genes > 10, 12, 
                                           ifelse(num_genes > 5, 14, 16)))))
  
  return(list(dynamic_space = dynamic_space, plot_height = plot_height, font_size = font_size))
}

create_heatmap_plot <- function(normalized_data, col_annotations, dynamic_space, plot_height) {
  heatmap_plot <- plot_ly(
    x = colnames(normalized_data),
    y = rownames(normalized_data),
    z = as.matrix(normalized_data),
    width = "2000",
    height = plot_height,
    type = "heatmap",
    colors = colorRamp(c("blue", "white", "red")),
    showscale = TRUE,
    colorbar = list(
      title = list(
        text = "Min-Max Normalized Gene Expression",
        font = list(size = 14, family = "Arial", color = "black"),
        side = "right"
      ),
      x = 1.1,  # Moves the color bar to the right
      tickfont = list(size = 16)
    ),
    hovertemplate = paste(
      "Sample: %{x}<br>",
      "Gene: %{y}<br>",
      "Normalized Expression: %{z:.2f}<extra></extra>"
    )
  ) %>%
    add_annotations(
      x = colnames(normalized_data),
      y = 1.001,
      text = dynamic_space,
      showarrow = FALSE,
      font = list(size = 8, color = as.character(col_annotations)),
      bgcolor = col_annotations,
      xref = "x",
      yref = "paper",
      xanchor = "center",
      yanchor = "bottom",
      align = "center"
    )
  
  return(heatmap_plot)
}

calculate_y_positions <- function(num_items, y_start = 1.02, y_end = 0.6) {
  if (num_items == 1) {
    y_positions <- y_start
  } else {
    total_height <- y_start - y_end
    y_gap <- total_height / (num_items - 1)
    y_positions <- y_start - (seq_len(num_items) - 1) * y_gap
  }
  return(y_positions)
}

add_legends_to_plot <- function(plot, EXP_groups, col_annotations, x_position = 1.2) {
  unique_groups <- unique(EXP_groups)
  unique_colors <- unique(col_annotations)
  num_items <- length(unique_groups)
  
  # Compute y positions dynamically
  y_positions <- calculate_y_positions(num_items)
  
  exp_annotations <- lapply(seq_along(unique_groups), function(i) {
    list(
      x = x_position,  # Adjust x position
      y = y_positions[i],
      xref = "paper", 
      yref = "paper",
      text = unique_groups[i],
      showarrow = FALSE,
      font = list(color = "#000000", size = 11),
      bgcolor = unique_colors[i],
      borderpad = 4,
      bordercolor = "#FFFFFF",
      borderwidth = 1,
      align = "center",
      width = 100,
      height = 20
    )
  })
  
  plot <- plot %>%
    layout(
      xaxis = list(tickangle = -45, tickfont = list(size = 13)),
      margin = list(l = 100, r = 400, b = 50, t = 80),
      xaxis = list(scaleanchor = "x", scaleratio = 1),
      yaxis = list(scaleanchor = "y", scaleratio = 1),
      legend = list(
        itemsizing = "constant",
        traceorder = "normal",
        font = list(family = "sans-serif", size = 12, color = "#000"),
        bgcolor = "#E2E2E2",
        bordercolor = "#FFFFFF",
        borderwidth = 2,
        orientation = "v"
      ),
      annotations = exp_annotations
    )
  
  return(plot)
}

kmeans_heatmap <- function(HM_expression_data_with_cluster_info, grps, Group_or_Sample) {
  # Prepare data and remove clustering columns
  HM_expression_data_ordered <- prepare_heatmap_data(
    HM_expression_data_with_cluster_info, grps, Group_or_Sample, exclude_cols = c("Cluster", "ClusterColor")
  )
  
  # Extract cluster information
  cluster_info <- HM_expression_data_with_cluster_info %>% select(Cluster, ClusterColor)
  
  # Generate annotations
  annotations_list <- generate_heatmap_annotations(grps, Group_or_Sample)
  EXP_groups <- annotations_list$EXP_groups
  col_annotations <- annotations_list$col_annotations
  
  num_genes <- nrow(HM_expression_data_ordered)
  num_columns <- ncol(HM_expression_data_ordered)
  
  # Calculate dynamic parameters
  params <- calculate_heatmap_parameters(num_genes, num_columns)
  dynamic_space <- params$dynamic_space
  plot_height <- params$plot_height
  font_size <- params$font_size
  
  # Create the base heatmap
  p <- create_heatmap_plot(HM_expression_data_ordered, col_annotations, dynamic_space, plot_height)
  
  # Add cluster annotations
  p <- p %>%
    add_trace(
      y = seq_along(cluster_info$Cluster),
      x = rep(ncol(HM_expression_data_ordered) + 1, length(cluster_info$Cluster)),
      z = matrix(cluster_info$ClusterColor, nrow = length(cluster_info$Cluster), ncol = 1),
      type = "heatmap",
      colorscale = "Greys",
      showscale = FALSE
    ) %>%
    layout(
      annotations = lapply(seq_along(cluster_info$Cluster), function(i) {
        list(
          x = 1.05,
          y = i - 1,
          text = paste0("Cluster ", cluster_info$Cluster[i]),
          xref = "paper",
          yref = "y",
          showarrow = FALSE,
          font = list(size = font_size),
          bgcolor = cluster_info$ClusterColor[i]
        )
      })
    )
  
  # Add legends
  p <- add_legends_to_plot(p, EXP_groups, col_annotations, x_position = 1.25)
  
  return(p)
}

HierClust_HM <- function(clean_variable_genes, grps, Group_or_Sample) {
  # Prepare data
  clean_variable_genes_ordered <- prepare_heatmap_data(clean_variable_genes, grps, Group_or_Sample)
  
  # Perform hierarchical clustering
  hclust_rows <- hierarchical_clustering(clean_variable_genes_ordered)
  ordered_GE_data <- clean_variable_genes_ordered[hclust_rows$order, ]
  
  # Min-max normalization
  normalized_data_t <- apply(ordered_GE_data, 1, min_max_normalization, min_value = -1, max_value = 1)
  normalized_data <- t(as.data.frame(normalized_data_t))
  
  # Generate annotations
  annotations_list <- generate_heatmap_annotations(grps, Group_or_Sample)
  EXP_groups <- annotations_list$EXP_groups
  col_annotations <- annotations_list$col_annotations
  
  num_genes <- nrow(normalized_data)
  num_columns <- ncol(normalized_data)
  
  # Calculate dynamic parameters
  params <- calculate_heatmap_parameters(num_genes, num_columns)
  dynamic_space <- params$dynamic_space
  plot_height <- params$plot_height
  font_size <- params$font_size
  
  # Create the base heatmap
  heatmap_plot <- create_heatmap_plot(normalized_data, col_annotations, dynamic_space, plot_height)
  
  # Add legends
  heatmap_plot <- add_legends_to_plot(heatmap_plot, EXP_groups, col_annotations, x_position = 1.28)
  
  return(heatmap_plot)
}

no_clustering_heatmap <- function(clean_variable_genes, grps, Group_or_Sample) {
  # Prepare data
  clean_variable_genes_ordered <- prepare_heatmap_data(clean_variable_genes, grps, Group_or_Sample)
  
  # Min-max normalization
  normalized_data_t <- apply(clean_variable_genes_ordered, 1, min_max_normalization, min_value = -1, max_value = 1)
  normalized_data <- t(as.data.frame(normalized_data_t))
  
  # Generate annotations
  annotations_list <- generate_heatmap_annotations(grps, Group_or_Sample)
  EXP_groups <- annotations_list$EXP_groups
  col_annotations <- annotations_list$col_annotations
  
  num_genes <- nrow(normalized_data)
  num_columns <- ncol(normalized_data)
  
  # Calculate dynamic parameters
  params <- calculate_heatmap_parameters(num_genes, num_columns)
  dynamic_space <- params$dynamic_space
  plot_height <- params$plot_height
  font_size <- params$font_size
  
  # Create the base heatmap
  heatmap_plot <- create_heatmap_plot(normalized_data, col_annotations, dynamic_space, plot_height)
  
  # Add legends
  heatmap_plot <- add_legends_to_plot(heatmap_plot, EXP_groups, col_annotations, x_position = 1.2)
  
  return(heatmap_plot)
}