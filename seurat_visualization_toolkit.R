# seurat_visualization_toolkit.R
# A Seurat Visualization Toolkit
# Author: Hanzhang
# Version: 2.0.0

# Data Flow Design Principles:
# 1. Data -> Extract coordinates -> Build layers -> Apply styles -> Output
# 2. Each function does one thing only
# 3. Use data structures rather than parameters to control behavior

library(Seurat)
library(ggplot2)
library(dplyr)

# ========== Core Data Structures ==========

#' Create visualization configuration
new_viz_config <- function(colors = NULL, point_size = 0.8, theme = "clean") {
  structure(
    list(
      colors = colors,
      point_size = point_size,
      theme = theme
    ),
    class = "viz_config"
  )
}

#' Extract plotting data from Seurat object
extract_plot_data <- function(seurat_obj, reduction, group_by) {
  # Invariant: ensure data integrity
  stopifnot(
    inherits(seurat_obj, "Seurat"),
    reduction %in% names(seurat_obj@reductions),
    group_by %in% colnames(seurat_obj@meta.data)
  )
  
  coords <- Embeddings(seurat_obj, reduction = reduction)
  meta <- seurat_obj@meta.data[, group_by, drop = FALSE]
  
  # Return standardized data structure
  cbind(
    data.frame(
      x = coords[, 1],
      y = coords[, 2]
    ),
    group = factor(meta[[group_by]])
  )
}

# ========== Color System ==========

# Predefined color palettes
PALETTES <- list(
  default = c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#9B59B6",
              "#1ABC9C", "#34495E", "#E67E22", "#95A5A6", "#D35400"),
  bright = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FECA57",
             "#F38BA8", "#A8E6CF", "#FFD93D", "#6BCF7F", "#4D96FF"),
  muted = c("#7F8C8D", "#95A5A6", "#BDC3C7", "#D5DBDB", "#E8DAEF",
            "#D6EAF8", "#D5F4E6", "#FCF3CF", "#FADBD8", "#EBDEF0")
)

#' Get color palette
get_colors <- function(n, palette = "default") {
  base_colors <- PALETTES[[palette]] %||% PALETTES$default
  
  if (n <= length(base_colors)) {
    base_colors[seq_len(n)]
  } else {
    colorRampPalette(base_colors)(n)
  }
}

# ========== Theme System ==========

#' Clean theme for general use
theme_clean <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", size = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

#' Publication-ready theme
theme_publication <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"),
      legend.key = element_blank()
    )
}

# Theme registry
THEMES <- list(
  clean = theme_clean,
  publication = theme_publication
)

#' Get theme by name
get_theme <- function(name = "clean") {
  theme_func <- THEMES[[name]] %||% THEMES$clean
  theme_func()
}

# ========== Core Visualization Functions ==========

#' Create dimensional reduction scatter plot
plot_reduction <- function(seurat_obj, 
                           reduction = "umap", 
                           group_by = "seurat_clusters",
                           config = new_viz_config()) {
  
  # Extract data
  plot_data <- extract_plot_data(seurat_obj, reduction, group_by)
  
  # Determine colors
  n_groups <- nlevels(plot_data$group)
  colors <- config$colors %||% get_colors(n_groups)
  
  # Build plot
  ggplot(plot_data, aes(x = x, y = y, color = group)) +
    geom_point(size = config$point_size, alpha = 0.8) +
    scale_color_manual(values = colors, name = group_by) +
    get_theme(config$theme) +
    labs(
      title = paste(toupper(reduction), "Plot"),
      x = paste(toupper(reduction), "1"),
      y = paste(toupper(reduction), "2")
    ) +
    coord_equal()
}

#' Add cluster labels to plot
add_cluster_labels <- function(plot, plot_data, repel = TRUE) {
  # Calculate label positions
  label_data <- plot_data %>%
    group_by(group) %>%
    summarise(x = median(x), y = median(y), .groups = 'drop')
  
  if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
    plot + ggrepel::geom_text_repel(
      data = label_data,
      aes(x = x, y = y, label = group),
      color = "black",
      fontface = "bold",
      size = 3.5,
      inherit.aes = FALSE
    )
  } else {
    plot + geom_text(
      data = label_data,
      aes(x = x, y = y, label = group),
      color = "black",
      fontface = "bold",
      size = 3.5,
      inherit.aes = FALSE
    )
  }
}

#' Create feature expression plots
plot_features <- function(seurat_obj, 
                          features, 
                          reduction = "umap",
                          config = new_viz_config()) {
  
  # Validate features exist
  available_features <- intersect(features, rownames(seurat_obj))
  if (length(available_features) == 0) {
    stop("No valid features found")
  }
  
  # Extract coordinates and expression data
  coords <- Embeddings(seurat_obj, reduction = reduction)
  expr_data <- FetchData(seurat_obj, vars = available_features)
  
  plot_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    expr_data,
    check.names = FALSE
  )
  
  # Create plot list
  plots <- map(available_features, function(feature) {
    ggplot(plot_data, aes(x = x, y = y)) +
      geom_point(aes(color = .data[[feature]]), 
                 size = config$point_size, alpha = 0.8) +
      scale_color_gradient2(
        low = "lightgrey", 
        mid = "yellow", 
        high = "red",
        midpoint = 0,
        name = "Expression"
      ) +
      get_theme(config$theme) +
      labs(
        title = feature,
        x = paste(toupper(reduction), "1"),
        y = paste(toupper(reduction), "2")
      ) +
      coord_equal()
  })
  
  names(plots) <- available_features
  plots
}

#' Create violin plots
plot_violin <- function(seurat_obj, 
                        features, 
                        group_by = "seurat_clusters",
                        config = new_viz_config()) {
  
  # Extract data
  plot_data <- FetchData(seurat_obj, vars = c(features, group_by))
  
  # Convert to long format
  plot_data_long <- tidyr::pivot_longer(
    plot_data, 
    cols = all_of(features),
    names_to = "feature", 
    values_to = "expression"
  )
  
  # Set colors
  n_groups <- length(unique(plot_data[[group_by]]))
  colors <- config$colors %||% get_colors(n_groups)
  
  ggplot(plot_data_long, aes(x = .data[[group_by]], y = expression)) +
    geom_violin(aes(fill = .data[[group_by]]), trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = colors, guide = "none") +
    get_theme(config$theme) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = group_by, y = "Expression") +
    facet_wrap(~ feature, scales = "free_y")
}

# ========== Convenience Interface ==========

#' Quick cluster plot creation
quick_cluster_plot <- function(seurat_obj, reduction = "umap", with_labels = TRUE) {
  plot_data <- extract_plot_data(seurat_obj, reduction, "seurat_clusters")
  
  p <- plot_reduction(seurat_obj, reduction = reduction)
  
  if (with_labels) {
    p <- add_cluster_labels(p, plot_data)
  }
  
  p
}

#' Quick feature plot creation
quick_feature_plot <- function(seurat_obj, features, reduction = "umap") {
  plots <- plot_features(seurat_obj, features, reduction)
  
  if (length(plots) == 1) {
    plots[[1]]
  } else if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(plots)
  } else {
    plots
  }
}

# ========== Utility Functions ==========

#' Save plot to file
save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  ggsave(filename, plot, width = width, height = height, dpi = dpi)
  message("Plot saved to: ", filename)
}

#' Get cell coordinates from reduction
get_coordinates <- function(seurat_obj, reduction = "umap") {
  coords <- Embeddings(seurat_obj, reduction = reduction)
  data.frame(
    cell = rownames(coords),
    x = coords[, 1],
    y = coords[, 2]
  )
}

# ========== Helper Functions ==========

#' Null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Simple map function
map <- function(.x, .f, ...) {
  lapply(.x, .f, ...)
}

# Usage examples
message("Seurat Visualization Toolkit v2.0 loaded!")
message("Quick start:")
message("  quick_cluster_plot(seurat_obj)")
message("  quick_feature_plot(seurat_obj, c('gene1', 'gene2'))")
