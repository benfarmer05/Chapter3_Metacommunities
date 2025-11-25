  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(sf)
  library(ggplot2)
  library(leaflet)
  library(cowplot)
  library(extrafont)
  library(tidyterra)
  library(ggnewscale)
  
  # For reading MATLAB v7.3 files
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("rhdf5")
  }
  library(rhdf5)
  
  ################################## PLOTTING OPTIONS ##################################
  
  # Font settings
  FONT_FAMILY <- "Georgia"
  
  # Susceptibility colors
  SUSC_COLORS <- c(
    "LS" = "#1E90FF",
    "MS" = "#FFD700",
    "HS" = "#FF1493"
  )
  
  # Plot bounds (single plot, no subplots)
  PLOT_BOUNDS <- list(
    xlim = c(-65.50911, -64.11830),
    ylim = c(18.04508, 18.83152)
  )
  
  # Figure dimensions
  FIG_WIDTH <- 10
  FIG_HEIGHT <- 5
  FIG_DPI <- 300
  
  # Text sizes
  TITLE_SIZE <- 12
  AXIS_TEXT_SIZE <- 7
  LABEL_TEXT_SIZE <- 2.5
  LEGEND_TEXT_SIZE <- 8
  LEGEND_TITLE_SIZE <- 9
  
  # Bubble plot settings
  # BUBBLE_SIZE_RANGE <- c(0.2, 2.5)  # min and max bubble sizes
  BUBBLE_SIZE_RANGE <- c(0.1, 2)  # min and max bubble sizes
  BUBBLE_ALPHA <- 0.7
  BUBBLE_STROKE <- 0.2  # outline thickness
  
  # Background colors
  PANEL_BG_COLOR <- "white"
  OCEAN_BG_COLOR <- "gray25"  # For areas with no data (outside depth range)
  
  # Island labels
  LABELS_ISLANDS <- data.frame(
    name = c("St. Thomas", "St. John", "Culebra", "Vieques", "Tortola", "Anegada"),
    x = c(-64.93, -64.75, -65.29, -65.40414, -64.606062, -64.35359),
    y = c(18.35, 18.337, 18.32, 18.13674, 18.438058, 18.736)
  )
  
  ################################## Load centroid data and grid ##################################
  
  cat("Loading centroid data...\n")
  centroids_data <- read.csv(here("data", "centroids_vertices_FINALFORCMS.csv"))
  
  cat("Initial reef sites:", nrow(centroids_data), "\n")
  
  cat("Creating grid polygons from vertices...\n")
  grid_polys <- lapply(1:nrow(centroids_data), function(i) {
    coords <- matrix(c(
      centroids_data$v1_lon[i], centroids_data$v1_lat[i],
      centroids_data$v2_lon[i], centroids_data$v2_lat[i],
      centroids_data$v3_lon[i], centroids_data$v3_lat[i],
      centroids_data$v4_lon[i], centroids_data$v4_lat[i],
      centroids_data$v1_lon[i], centroids_data$v1_lat[i]
    ), ncol = 2, byrow = TRUE)
    st_polygon(list(coords))
  })
  
  grid_sf <- st_sf(
    unique_ID = centroids_data$unique_ID,
    mean_coral_cover = centroids_data$mean_coral_cover,
    low_coral_cover = centroids_data$low_coral_cover,
    moderate_coral_cover = centroids_data$moderate_coral_cover,
    high_coral_cover = centroids_data$high_coral_cover,
    centroid_lon = centroids_data$centroid_lon,
    centroid_lat = centroids_data$centroid_lat,
    geometry = st_sfc(grid_polys, crs = 4326)
  )
  
  cat("Grid polygons created\n\n")
  
  ################################## Load MATLAB connectivity results ##################################
  
  cat("Loading connectivity matrices from MATLAB...\n")
  
  matlab_cache_file <- here("temp", "conn_structs_FINAL_20181201_to_20191231_notraj_2018bf.mat")
  
  if (!file.exists(matlab_cache_file)) {
    stop("MATLAB cache file not found: ", matlab_cache_file, 
         "\nPlease run the MATLAB script with RECALCULATE_CONNECTIVITY = true first.")
  }
  
  h5_contents <- h5ls(matlab_cache_file)
  original_indices <- as.vector(h5read(matlab_cache_file, "/final_indices"))
  
  cat("Connectivity matrices use:", length(original_indices), "reefs\n")
  cat("Original indices range:", range(original_indices), "\n")
  cat("Full reef dataset has:", nrow(centroids_data), "reefs\n\n")
  
  ################################## Identify non-source reefs ##################################
  
  cat("Identifying reefs that never produced trajectories...\n")
  
  n_filtered_reefs <- length(original_indices)
  filtered_sources_ever <- rep(FALSE, n_filtered_reefs)
  
  refs_contents <- h5_contents[h5_contents$group == "/#refs#" & h5_contents$otype == "H5I_GROUP", ]
  sparse_matrix_names <- refs_contents$name
  
  cat("Found", length(sparse_matrix_names), "sparse matrices in /#refs#\n")
  
  for (i in seq_along(sparse_matrix_names)) {
    matrix_name <- sparse_matrix_names[i]
    matrix_path <- paste0("/#refs#/", matrix_name)
    
    tryCatch({
      ir <- h5read(matlab_cache_file, paste0(matrix_path, "/ir"))
      if (length(ir) > 0) {
        row_indices <- ir + 1
        filtered_sources_ever[unique(row_indices)] <- TRUE
      }
    }, error = function(e) {})
    
    if (i %% 10 == 0) cat("  Processed", i, "of", length(sparse_matrix_names), "matrices\n")
  }
  
  cat("  Processed all", length(sparse_matrix_names), "matrices\n")
  
  n_reefs <- nrow(centroids_data)
  all_sources_ever <- rep(FALSE, n_reefs)
  reefs_in_filtered_set <- original_indices
  all_sources_ever[reefs_in_filtered_set] <- filtered_sources_ever
  
  cat("\nDiagnostic info:\n")
  cat("  Total reefs in full dataset:", n_reefs, "\n")
  cat("  Reefs in connectivity matrices:", n_filtered_reefs, "\n")
  cat("  Of those, reefs with trajectories:", sum(filtered_sources_ever), "\n")
  cat("  Overall, reefs marked as sources:", sum(all_sources_ever), "\n")
  
  inactive_sources <- which(!all_sources_ever)
  active_sources <- which(all_sources_ever)
  
  cat("\n=== TRAJECTORY ANALYSIS RESULTS ===\n")
  cat("Reefs that NEVER produced trajectories:", length(inactive_sources), 
      sprintf("(%.1f%%)\n", 100 * length(inactive_sources) / n_reefs))
  cat("Reefs with trajectories:", length(active_sources), 
      sprintf("(%.1f%%)\n", 100 * length(active_sources) / n_reefs))
  cat("====================================\n\n")
  
  ################################## Filter data ##################################
  
  cat("\nFiltering data to active sources only...\n")
  
  centroids_filtered <- centroids_data[active_sources, ]
  grid_filtered_sf <- grid_sf[active_sources, ]
  grid_inactive_sf <- grid_sf[inactive_sources, ]
  
  cat("Filtered reef sites:", nrow(centroids_filtered), "\n")
  cat("Reduction:", round((1 - nrow(centroids_filtered)/nrow(centroids_data)) * 100, 1), "%\n\n")
  
  ################################## Create Dominant Susceptibility ##################################
  
  cat("Creating dominant susceptibility classification...\n")
  
  grid_filtered_sf$dominant_susc <- NA_character_
  
  ls_mask <- grid_filtered_sf$low_coral_cover >= grid_filtered_sf$moderate_coral_cover & 
    grid_filtered_sf$low_coral_cover >= grid_filtered_sf$high_coral_cover
  grid_filtered_sf$dominant_susc[ls_mask] <- "LS"
  
  ms_mask <- grid_filtered_sf$moderate_coral_cover > grid_filtered_sf$low_coral_cover & 
    grid_filtered_sf$moderate_coral_cover >= grid_filtered_sf$high_coral_cover
  grid_filtered_sf$dominant_susc[ms_mask] <- "MS"
  
  hs_mask <- grid_filtered_sf$high_coral_cover > grid_filtered_sf$low_coral_cover & 
    grid_filtered_sf$high_coral_cover > grid_filtered_sf$moderate_coral_cover
  grid_filtered_sf$dominant_susc[hs_mask] <- "HS"
  
  grid_filtered_sf$dominant_susc <- factor(grid_filtered_sf$dominant_susc, levels = c("LS", "MS", "HS"))
  
  cat("Dominant susceptibility calculated\n")
  cat("  LS dominant:", sum(grid_filtered_sf$dominant_susc == "LS", na.rm = TRUE), "sites\n")
  cat("  MS dominant:", sum(grid_filtered_sf$dominant_susc == "MS", na.rm = TRUE), "sites\n")
  cat("  HS dominant:", sum(grid_filtered_sf$dominant_susc == "HS", na.rm = TRUE), "sites\n\n")
  
  ################################## Spatial objects for mapping ##################################
  
  cat("Creating spatial objects for mapping...\n")
  grid_sf$has_trajectories <- all_sources_ever
  cat("Spatial objects created\n\n")
  
  ################################## Load OSM land data ##################################
  
  cat("Loading land polygons...\n")
  if (file.exists(here("output", "osm_land_vect.rds"))) {
    land_vect <- readRDS(here("output", "osm_land_vect.rds"))
    land_sf <- st_as_sf(land_vect)
    cat("Land polygons loaded\n\n")
  } else {
    land_sf <- NULL
    cat("No land polygons found - maps will not show land\n\n")
  }
  
  ################################## Load bathymetry and create depth mask ##################################
  
  cat("Loading bathymetry data for depth masking...\n")
  
  spatial_metadata <- readRDS(here('data/spatial_metadata.rds'))
  bathy_final <- readRDS(here('data/bathy_final.rds'))
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  
  # Create depth mask for 0-60m (bathy values are negative, so 0 to -60)
  depth_mask <- bathy_final >= -60 & bathy_final <= 0
  
  # Create a raster that is gray (value = 1) within depth range, NA outside
  depth_zone_raster <- depth_mask
  depth_zone_raster[depth_mask == 0] <- NA
  depth_zone_raster[depth_mask == 1] <- 1
  
  # Project to lat/lon for plotting
  depth_zone_latlon <- project(depth_zone_raster, "EPSG:4326", method = "near")
  
  cat("Depth mask created (0-60m zone)\n\n")
  
  ################################## Create leaflet maps ##################################
  
  # cat("Creating interactive leaflet maps...\n\n")
  # 
  # pal_total <- colorNumeric(palette = "viridis", domain = grid_filtered_sf$mean_coral_cover * 100)
  # pal_low <- colorNumeric(palette = "plasma", domain = grid_filtered_sf$low_coral_cover * 100)
  # pal_moderate <- colorNumeric(palette = "inferno", domain = grid_filtered_sf$moderate_coral_cover * 100)
  # pal_high <- colorNumeric(palette = "magma", domain = grid_filtered_sf$high_coral_cover * 100)
  # 
  # cat("Map 1: Total coral cover\n")
  # map1 <- leaflet() %>%
  #   addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  #   addTiles(group = "OpenStreetMap") %>%
  #   addPolygons(data = grid_filtered_sf,
  #               fillColor = ~pal_total(mean_coral_cover * 100),
  #               fillOpacity = 0.6, color = "black", weight = 1,
  #               popup = ~paste0("<b>ID:</b> ", unique_ID, 
  #                               "<br><b>Total Cover:</b> ", round(mean_coral_cover * 100, 2), "%")) %>%
  #   addLegend("bottomright", pal = pal_total, 
  #             values = grid_filtered_sf$mean_coral_cover * 100,
  #             title = "Total Coral<br>Cover (%)", opacity = 1) %>%
  #   addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
  #                    options = layersControlOptions(collapsed = FALSE))
  # print(map1)
  # 
  # cat("Map 2: Low susceptibility cover\n")
  # map2 <- leaflet() %>%
  #   addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  #   addTiles(group = "OpenStreetMap") %>%
  #   addPolygons(data = grid_filtered_sf,
  #               fillColor = ~pal_low(low_coral_cover * 100),
  #               fillOpacity = 0.6, color = "black", weight = 1,
  #               popup = ~paste0("<b>ID:</b> ", unique_ID, 
  #                               "<br><b>Low Susceptibility:</b> ", round(low_coral_cover * 100, 2), "%")) %>%
  #   addLegend("bottomright", pal = pal_low, 
  #             values = grid_filtered_sf$low_coral_cover * 100,
  #             title = "Low Susceptibility<br>Cover (%)", opacity = 1) %>%
  #   addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
  #                    options = layersControlOptions(collapsed = FALSE))
  # print(map2)
  # 
  # cat("Map 3: Moderate susceptibility cover\n")
  # map3 <- leaflet() %>%
  #   addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  #   addTiles(group = "OpenStreetMap") %>%
  #   addPolygons(data = grid_filtered_sf,
  #               fillColor = ~pal_moderate(moderate_coral_cover * 100),
  #               fillOpacity = 0.6, color = "black", weight = 1,
  #               popup = ~paste0("<b>ID:</b> ", unique_ID, 
  #                               "<br><b>Moderate Susceptibility:</b> ", round(moderate_coral_cover * 100, 2), "%")) %>%
  #   addLegend("bottomright", pal = pal_moderate, 
  #             values = grid_filtered_sf$moderate_coral_cover * 100,
  #             title = "Moderate Susceptibility<br>Cover (%)", opacity = 1) %>%
  #   addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
  #                    options = layersControlOptions(collapsed = FALSE))
  # print(map3)
  # 
  # cat("Map 4: High susceptibility cover\n")
  # map4 <- leaflet() %>%
  #   addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  #   addTiles(group = "OpenStreetMap") %>%
  #   addPolygons(data = grid_filtered_sf,
  #               fillColor = ~pal_high(high_coral_cover * 100),
  #               fillOpacity = 0.6, color = "black", weight = 1,
  #               popup = ~paste0("<b>ID:</b> ", unique_ID, 
  #                               "<br><b>High Susceptibility:</b> ", round(high_coral_cover * 100, 2), "%")) %>%
  #   addLegend("bottomright", pal = pal_high, 
  #             values = grid_filtered_sf$high_coral_cover * 100,
  #             title = "High Susceptibility<br>Cover (%)", opacity = 1) %>%
  #   addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
  #                    options = layersControlOptions(collapsed = FALSE))
  # print(map4)
  # 
  # cat("Map 5: Active vs. Inactive source reefs\n")
  # map5 <- leaflet() %>%
  #   addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  #   addTiles(group = "OpenStreetMap") %>%
  #   addPolygons(data = grid_inactive_sf,
  #               fillColor = "red", fillOpacity = 0.6, color = "darkred", weight = 1,
  #               popup = ~paste0("<b>ID:</b> ", unique_ID, "<br><b>Status:</b> No trajectories")) %>%
  #   addPolygons(data = grid_filtered_sf,
  #               fillColor = ~pal_total(mean_coral_cover * 100),
  #               fillOpacity = 0.6, color = "black", weight = 1,
  #               popup = ~paste0("<b>ID:</b> ", unique_ID, 
  #                               "<br><b>Total Cover:</b> ", round(mean_coral_cover * 100, 2), "%",
  #                               "<br><b>Status:</b> Active source")) %>%
  #   addLegend("bottomright", colors = c("red", "blue"),
  #             labels = c(sprintf("No trajectories (n=%d)", length(inactive_sources)),
  #                        sprintf("Active sources (n=%d)", nrow(grid_filtered_sf))),
  #             title = "Reef Status", opacity = 1) %>%
  #   addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
  #                    options = layersControlOptions(collapsed = FALSE))
  # print(map5)
  # 
  # cat("Map 6: Dominant susceptibility (leaflet)\n")
  # pal_susc <- colorFactor(palette = unname(SUSC_COLORS), levels = c("LS", "MS", "HS"))
  # 
  # map6 <- leaflet() %>%
  #   addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  #   addTiles(group = "OpenStreetMap") %>%
  #   addPolygons(data = grid_filtered_sf,
  #               fillColor = ~pal_susc(dominant_susc),
  #               fillOpacity = 0.7, color = "black", weight = 0.5,
  #               popup = ~paste0("<b>ID:</b> ", unique_ID,
  #                               "<br><b>Dominant:</b> ", dominant_susc,
  #                               "<br><b>LS Cover:</b> ", round(low_coral_cover * 100, 2), "%",
  #                               "<br><b>MS Cover:</b> ", round(moderate_coral_cover * 100, 2), "%",
  #                               "<br><b>HS Cover:</b> ", round(high_coral_cover * 100, 2), "%",
  #                               "<br><b>Total Cover:</b> ", round(mean_coral_cover * 100, 2), "%")) %>%
  #   addLegend("bottomright", pal = pal_susc, values = c("LS", "MS", "HS"),
  #             title = "Dominant<br>Susceptibility",
  #             labels = c("Low (LS)", "Moderate (MS)", "High (HS)"), opacity = 1) %>%
  #   addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
  #                    options = layersControlOptions(collapsed = FALSE))
  # print(map6)
  
  ################################## ggplot: Dominant Susceptibility Polygon Map ##################################
  
  cat("\nCreating ggplot dominant susceptibility map (polygon style)...\n")
  
  create_susc_polygon_plot <- function(bounds, data_sf, land_data = NULL, 
                                       depth_raster = NULL, labels_df = NULL, show_legend = TRUE) {
    
    # Crop depth raster to bounds
    buffer <- 0.1
    region_extent <- ext(bounds$xlim[1]-buffer, bounds$xlim[2]+buffer, 
                         bounds$ylim[1]-buffer, bounds$ylim[2]+buffer)
    
    p <- ggplot()
    
    # Add depth zone background first
    if (!is.null(depth_raster)) {
      depth_cropped <- crop(depth_raster, region_extent)
      p <- p + geom_spatraster(data = depth_cropped) +
        scale_fill_gradient(low = OCEAN_BG_COLOR, high = OCEAN_BG_COLOR, na.value = "transparent", guide = "none")
    }
    
    # Add susceptibility polygons on top
    p <- p + 
      ggnewscale::new_scale_fill() +
      geom_sf(data = data_sf, aes(fill = dominant_susc), color = NA) +
      scale_fill_manual(values = SUSC_COLORS,
                        labels = c("Low", "Moderate", "High"),
                        name = "Dominant\nSusceptibility",
                        drop = FALSE)
    
    if (!is.null(land_data)) {
      p <- p + geom_sf(data = land_data, fill = "gray90", color = "black", linewidth = 0.3)
    }
    
    # Add island labels
    if (!is.null(labels_df) && nrow(labels_df) > 0) {
      p <- p + geom_text(data = labels_df, aes(x = x, y = y, label = name),
                         size = LABEL_TEXT_SIZE, family = FONT_FAMILY, color = "black")
    }
    
    p <- p +
      coord_sf(xlim = bounds$xlim, ylim = bounds$ylim, expand = FALSE) +
      theme_classic(base_family = FONT_FAMILY) +
      theme(
        legend.position = if(show_legend) "right" else "none",
        legend.text = element_text(size = LEGEND_TEXT_SIZE),
        legend.title = element_text(size = LEGEND_TITLE_SIZE),
        axis.title = element_blank(),
        axis.text = element_text(size = AXIS_TEXT_SIZE),
        plot.margin = margin(5, 5, 5, 5, 'pt'),
        panel.background = element_rect(fill = PANEL_BG_COLOR, color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    return(p)
  }
  
  gg_susc_poly <- create_susc_polygon_plot(PLOT_BOUNDS, grid_filtered_sf, land_sf, 
                                           depth_zone_latlon, LABELS_ISLANDS, TRUE)
  print(gg_susc_poly)
  
  ggsave(filename = here("output", "output_figures_tables", "fig_dominant_susceptibility_polygons.png"),
         plot = gg_susc_poly, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI, bg = "white")
  cat("Saved: output_figures_tables/fig_dominant_susceptibility_polygons.png\n")
  
  ################################## ggplot: Bubble Plot with Coral Cover ##################################
  
  cat("\nCreating ggplot bubble plot (circles scaled by coral cover)...\n")
  
  # Create dataframe with centroids for bubble plot
  bubble_df <- data.frame(
    lon = grid_filtered_sf$centroid_lon,
    lat = grid_filtered_sf$centroid_lat,
    total_cover = grid_filtered_sf$mean_coral_cover,
    dominant_susc = grid_filtered_sf$dominant_susc
  )
  
  # Remove rows with zero cover
  bubble_df_nonzero <- bubble_df[bubble_df$total_cover > 0, ]
  
  # Calculate max cover for scaling (as percentage)
  max_cover <- max(bubble_df_nonzero$total_cover, na.rm = TRUE)
  cat("  Max coral cover:", round(max_cover * 100, 1), "%\n")
  
  create_susc_bubble_plot <- function(bounds, data_df, land_data = NULL, 
                                      depth_raster = NULL, labels_df = NULL, show_legend = TRUE) {
    
    # Crop depth raster to bounds
    buffer <- 0.1
    region_extent <- ext(bounds$xlim[1]-buffer, bounds$xlim[2]+buffer, 
                         bounds$ylim[1]-buffer, bounds$ylim[2]+buffer)
    
    # Filter to bounds
    data_region <- data_df[data_df$lon >= bounds$xlim[1] & data_df$lon <= bounds$xlim[2] &
                             data_df$lat >= bounds$ylim[1] & data_df$lat <= bounds$ylim[2], ]
    
    p <- ggplot()
    
    # Add depth zone background first
    if (!is.null(depth_raster)) {
      depth_cropped <- crop(depth_raster, region_extent)
      p <- p + geom_spatraster(data = depth_cropped) +
        scale_fill_gradient(low = OCEAN_BG_COLOR, high = OCEAN_BG_COLOR, na.value = "transparent", guide = "none")
    }
    
    if (!is.null(land_data)) {
      p <- p + geom_sf(data = land_data, fill = "gray90", color = "black", linewidth = 0.3)
    }
    
    # Bubbles with outline only (no fill)
    p <- p +
      geom_point(data = data_region, 
                 aes(x = lon, y = lat, size = total_cover, color = dominant_susc),
                 shape = 21,  # Circle with outline
                 fill = NA,   # No fill
                 stroke = BUBBLE_STROKE,
                 alpha = BUBBLE_ALPHA) +
      scale_color_manual(values = SUSC_COLORS,
                         labels = c("Low", "Moderate", "High"),
                         name = "Susceptibility",
                         drop = FALSE) +
      scale_size_area(name = "Coral Cover",
                      max_size = BUBBLE_SIZE_RANGE[2],
                      limits = c(0, 0.5),
                      guide = "none") +
      guides(color = guide_legend(order = 1, override.aes = list(size = 4, stroke = 1.5)))
    
    # Add island labels
    if (!is.null(labels_df) && nrow(labels_df) > 0) {
      p <- p + geom_text(data = labels_df, aes(x = x, y = y, label = name),
                         size = LABEL_TEXT_SIZE, family = FONT_FAMILY, color = "black")
    }
    
    p <- p +
      coord_sf(xlim = bounds$xlim, ylim = bounds$ylim, expand = FALSE) +
      theme_classic(base_family = FONT_FAMILY) +
      theme(
        legend.position = if(show_legend) "right" else "none",
        legend.text = element_text(size = LEGEND_TEXT_SIZE),
        legend.title = element_text(size = LEGEND_TITLE_SIZE),
        axis.title = element_blank(),
        axis.text = element_text(size = AXIS_TEXT_SIZE),
        plot.margin = margin(5, 5, 5, 5, 'pt'),
        panel.background = element_rect(fill = PANEL_BG_COLOR, color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    return(p)
  }
  
  gg_susc_bubble <- create_susc_bubble_plot(PLOT_BOUNDS, bubble_df_nonzero, land_sf, 
                                            depth_zone_latlon, LABELS_ISLANDS, TRUE)
  print(gg_susc_bubble)
  
  ggsave(filename = here("output", "output_figures_tables", "fig_dominant_susceptibility_bubbles.png"),
         plot = gg_susc_bubble, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI, bg = "white")
  cat("Saved: output_figures_tables/fig_dominant_susceptibility_bubbles.png\n")
  
  ################################## Summary Statistics ##################################
  
  cat("\n=== SUMMARY STATISTICS ===\n")
  cat("Original reefs:", nrow(centroids_data), "\n")
  cat("Filtered reefs (active sources):", nrow(centroids_filtered), "\n")
  cat("Removed reefs (no trajectories):", length(inactive_sources), "\n")
  cat("Reduction: ", round((1 - nrow(centroids_filtered)/nrow(centroids_data)) * 100, 1), "%\n\n")
  
  cat("Mean coral cover (filtered reefs):\n")
  cat("  Total:", round(mean(centroids_filtered$mean_coral_cover) * 100, 2), "%\n")
  cat("  Low susceptibility:", round(mean(centroids_filtered$low_coral_cover) * 100, 2), "%\n")
  cat("  Moderate susceptibility:", round(mean(centroids_filtered$moderate_coral_cover) * 100, 2), "%\n")
  cat("  High susceptibility:", round(mean(centroids_filtered$high_coral_cover) * 100, 2), "%\n")
  
  cat("\n=== SUSCEPTIBILITY SUMMARY ===\n")
  
  susc_counts <- table(grid_filtered_sf$dominant_susc)
  cat("Sites by dominant susceptibility group:\n")
  for (g in c("LS", "MS", "HS")) {
    n <- ifelse(g %in% names(susc_counts), susc_counts[g], 0)
    pct <- 100 * n / nrow(grid_filtered_sf)
    cat(sprintf("  %s: %d sites (%.1f%%)\n", g, n, pct))
  }
  
  cat("\nTotal coral cover by susceptibility group:\n")
  total_ls <- sum(grid_filtered_sf$low_coral_cover)
  total_ms <- sum(grid_filtered_sf$moderate_coral_cover)
  total_hs <- sum(grid_filtered_sf$high_coral_cover)
  total_all <- total_ls + total_ms + total_hs
  
  cat(sprintf("  Low susceptibility:      %.2f (%.1f%% of total)\n", total_ls, 100 * total_ls / total_all))
  cat(sprintf("  Moderate susceptibility: %.2f (%.1f%% of total)\n", total_ms, 100 * total_ms / total_all))
  cat(sprintf("  High susceptibility:     %.2f (%.1f%% of total)\n", total_hs, 100 * total_hs / total_all))
  
  cat("==============================\n\n")
  
  cat("\n=== COMPLETE ===\n")