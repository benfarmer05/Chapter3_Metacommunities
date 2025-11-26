  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(rhdf5)
  library(terra)
  library(sf)
  library(tidyverse)
  library(tidyterra)
  library(patchwork)
  
  library(viridis)
  
  ################################## Settings ##################################
  
  # Plot styling parameters
  titlesize <- 10
  textsize <- 9
  axis_textsize <- 6
  legend_textsize <- 6
  legend_barwidth <- 0.3
  legend_barheight <- 0.5
  
  # Figure dimensions
  fig_width <- 6
  fig_height <- 3
  
  # Color palettes (viridis options: viridis, magma, plasma, inferno, cividis, mako, rocket, turbo)
  palette_I_abs <- "inferno"
  palette_I_rel <- "turbo"
  palette_R_abs <- "inferno"
  palette_R_rel <- "turbo"
  
  # Extent for plots
  obs_extent <- ext(-65.50911, -64.57022, 18.04508, 18.64853)
  
  # Disease threshold for "observable" removed cover
  removed_cover_thresh <- 0.001
  
  # Custom colormap matching MATLAB stacked colormap
  
  # Using viridis plasma instead - simpler and cleaner
  
  ################################## Load Data ##################################
  
  cat("Loading MATLAB workspace...\n")
  matlab_file <- here("output", "seascape_SIR", "SCENARIO1_ET0_I0_0p0000010455_FS1.0_FP1.50_better_removal", "seascape_SIR_workspace.mat")
  
  # Extract data
  locations <- h5read(matlab_file, "/outputs/sites/locations")
  if (ncol(locations) > nrow(locations)) locations <- t(locations)
  
  unique_IDs <- as.vector(h5read(matlab_file, "/outputs/sites/unique_IDs"))
  N_site <- as.vector(h5read(matlab_file, "/outputs/sites/N_site"))
  
  S_total <- h5read(matlab_file, "/outputs/totals/S_total")
  I_total <- h5read(matlab_file, "/outputs/totals/I_total")
  R_total <- h5read(matlab_file, "/outputs/totals/R_total")
  
  ref_date <- as.Date("2018-12-01")
  tspan <- as.vector(h5read(matlab_file, "/outputs/metadata/tspan"))
  sim_dates <- ref_date + (tspan[1]:tspan[2]) - 1
  
  num_days <- nrow(I_total)
  final_date <- max(sim_dates)
  
  cat("Loaded", ncol(I_total), "sites,", num_days, "days\n")
  cat("Final date:", as.character(final_date), "\n\n")
  
  ################################## Load Grid Polygons ##################################
  
  cat("Loading grid polygons...\n")
  centroids_data <- read.csv(here("data", "centroids_vertices_FINALFORCMS.csv"))
  land_vect <- readRDS(here("output", "osm_land_vect.rds"))
  
  # Only include grid cells that are in the model output (match by unique_ID)
  grid_data <- centroids_data[match(unique_IDs, centroids_data$unique_ID), ]
  grid_data$I_absolute <- I_total[num_days, ]
  grid_data$S_final <- S_total[num_days, ]
  grid_data$N_site <- N_site
  
  # Calculate values matching MATLAB approach
  grid_data$I_relative <- grid_data$I_absolute / grid_data$N_site
  grid_data$R_absolute <- grid_data$N_site - grid_data$S_final
  grid_data$R_relative <- (grid_data$N_site - grid_data$S_final) / grid_data$N_site
  
  # Create polygons
  grid_polys <- lapply(1:nrow(grid_data), function(i) {
    if (!is.na(grid_data$v1_lon[i])) {
      coords <- matrix(c(
        grid_data$v1_lon[i], grid_data$v1_lat[i],
        grid_data$v2_lon[i], grid_data$v2_lat[i],
        grid_data$v3_lon[i], grid_data$v3_lat[i],
        grid_data$v4_lon[i], grid_data$v4_lat[i],
        grid_data$v1_lon[i], grid_data$v1_lat[i]
      ), ncol = 2, byrow = TRUE)
      return(st_polygon(list(coords)))
    } else {
      return(NULL)
    }
  })
  
  valid_idx <- !sapply(grid_polys, is.null)
  grid_sf <- st_sf(
    unique_ID = grid_data$unique_ID[valid_idx],
    I_absolute = grid_data$I_absolute[valid_idx],
    I_relative = grid_data$I_relative[valid_idx],
    R_absolute = grid_data$R_absolute[valid_idx],
    R_relative = grid_data$R_relative[valid_idx],
    geometry = st_sfc(grid_polys[valid_idx], crs = 4326)
  )
  
  # Replace Inf/NaN with 0
  grid_sf$I_relative[!is.finite(grid_sf$I_relative)] <- 0
  grid_sf$R_relative[!is.finite(grid_sf$R_relative)] <- 0
  
  cat("Grid cells in model:", nrow(grid_sf), "\n")
  
  # Calculate max values for color scale limits (as %, rounded to 2 decimals)
  max_I_abs <- max(grid_sf$I_absolute, na.rm = TRUE) * 100
  max_I_rel <- max(grid_sf$I_relative, na.rm = TRUE) * 100
  max_R_abs <- max(grid_sf$R_absolute, na.rm = TRUE) * 100
  max_R_rel <- max(grid_sf$R_relative, na.rm = TRUE) * 100
  
  # Round up to 2 decimal places
  round_up_2 <- function(x) ceiling(x * 100) / 100
  max_I_abs_rounded <- round_up_2(max_I_abs)
  max_I_rel_rounded <- round_up_2(max_I_rel)
  max_R_abs_rounded <- round_up_2(max_R_abs)
  max_R_rel_rounded <- round_up_2(max_R_rel)
  
  cat("Max I absolute:", round(max_I_abs, 4), "% -> rounded to", max_I_abs_rounded, "%\n")
  cat("Max I relative:", round(max_I_rel, 4), "% -> rounded to", max_I_rel_rounded, "%\n")
  cat("Max R absolute:", round(max_R_abs, 4), "% -> rounded to", max_R_abs_rounded, "%\n")
  cat("Max R relative:", round(max_R_rel, 4), "% -> rounded to", max_R_rel_rounded, "%\n")
  
  ################################## Create Plots ##################################
  
  cat("Creating final state plots...\n")
  
  # Common theme
  theme_map <- theme_classic(base_family = "Georgia") +
    theme(
      axis.title = element_text(size = textsize, color = 'black'),
      axis.text = element_text(size = axis_textsize, color = 'black'),
      axis.ticks = element_line(color = "black"),
      legend.position = "right",
      legend.text = element_text(size = legend_textsize),
      legend.title = element_text(size = legend_textsize),
      legend.key.width = unit(legend_barwidth, "cm"),
      legend.key.height = unit(legend_barheight, "cm"),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  # Panel A: Absolute infected cover
  p_I_abs <- ggplot() +
    geom_sf(data = grid_sf, aes(fill = I_absolute * 100), color = NA) +
    scale_fill_viridis(name = "Abs. (%)", option = palette_I_abs, 
                       na.value = "transparent", limits = c(0, max_I_abs_rounded),
                       labels = function(x) sprintf("%.2f", x)) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "gray40", linewidth = 0.2) +
    coord_sf(xlim = c(obs_extent[1], obs_extent[2]), 
             ylim = c(obs_extent[3], obs_extent[4]), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    annotate("text", x = Inf, y = Inf, label = "A", hjust = 1.5, vjust = 1.5,
             size = titlesize / .pt, family = "Georgia") +
    theme_map
  
  # Panel B: Relative infected cover
  p_I_rel <- ggplot() +
    geom_sf(data = grid_sf, aes(fill = I_relative * 100), color = NA) +
    scale_fill_viridis(name = "Rel. (%)", option = palette_I_rel, 
                       na.value = "transparent", limits = c(0, max_I_rel_rounded),
                       labels = function(x) sprintf("%.2f", x)) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "gray40", linewidth = 0.2) +
    coord_sf(xlim = c(obs_extent[1], obs_extent[2]), 
             ylim = c(obs_extent[3], obs_extent[4]), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    annotate("text", x = Inf, y = Inf, label = "B", hjust = 1.5, vjust = 1.5,
             size = titlesize / .pt, family = "Georgia") +
    theme_map
  
  # Panel C: Absolute removed cover
  p_R_abs <- ggplot() +
    geom_sf(data = grid_sf, aes(fill = R_absolute * 100), color = NA) +
    scale_fill_viridis(name = "Abs. (%)", option = palette_R_abs, 
                       na.value = "transparent", limits = c(0, max_R_abs_rounded),
                       labels = function(x) sprintf("%.2f", x)) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "gray40", linewidth = 0.2) +
    coord_sf(xlim = c(obs_extent[1], obs_extent[2]), 
             ylim = c(obs_extent[3], obs_extent[4]), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    annotate("text", x = Inf, y = Inf, label = "C", hjust = 1.5, vjust = 1.5,
             size = titlesize / .pt, family = "Georgia") +
    theme_map
  
  # Panel D: Relative removed cover  
  p_R_rel <- ggplot() +
    geom_sf(data = grid_sf, aes(fill = R_relative * 100), color = NA) +
    scale_fill_viridis(name = "Rel. (%)", option = palette_R_rel, 
                       na.value = "transparent", limits = c(0, max_R_rel_rounded),
                       labels = function(x) sprintf("%.2f", x)) +
    geom_spatvector(data = land_vect, fill = "gray90", color = "gray40", linewidth = 0.2) +
    coord_sf(xlim = c(obs_extent[1], obs_extent[2]), 
             ylim = c(obs_extent[3], obs_extent[4]), expand = FALSE) +
    labs(x = NULL, y = NULL) +
    annotate("text", x = Inf, y = Inf, label = "D", hjust = 1.5, vjust = 1.5,
             size = titlesize / .pt, family = "Georgia") +
    theme_map
  
  # Combine plots
  final_state_plot <- (p_I_abs + p_I_rel) / (p_R_abs + p_R_rel)
  
  print(final_state_plot)
  
  # Save
  ggsave(
    filename = here("output", "output_figures_tables", "final_state_plot.png"),
    plot = final_state_plot,
    width = fig_width,
    height = fig_height,
    dpi = 300,
    bg = "white"
  )
  
  cat("\n✓ Final state plot saved!\n")