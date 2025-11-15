  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(terra)
  library(sf)
  library(ggplot2)
  library(leaflet)
  
  # For reading MATLAB v7.3 files
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("rhdf5")
  }
  library(rhdf5)
  
  ################################## Load centroid data and grid ##################################
  
  cat("Loading centroid data...\n")
  centroids_data <- read.csv(here("data", "centroids_vertices_FINALFORCMS.csv"))
  
  cat("Initial reef sites:", nrow(centroids_data), "\n")
  
  # Create grid polygons from vertices
  cat("Creating grid polygons from vertices...\n")
  grid_polys <- lapply(1:nrow(centroids_data), function(i) {
    coords <- matrix(c(
      centroids_data$v1_lon[i], centroids_data$v1_lat[i],
      centroids_data$v2_lon[i], centroids_data$v2_lat[i],
      centroids_data$v3_lon[i], centroids_data$v3_lat[i],
      centroids_data$v4_lon[i], centroids_data$v4_lat[i],
      centroids_data$v1_lon[i], centroids_data$v1_lat[i]  # Close the polygon
    ), ncol = 2, byrow = TRUE)
    st_polygon(list(coords))
  })
  
  grid_sf <- st_sf(
    unique_ID = centroids_data$unique_ID,
    mean_coral_cover = centroids_data$mean_coral_cover,
    low_coral_cover = centroids_data$low_coral_cover,
    moderate_coral_cover = centroids_data$moderate_coral_cover,
    high_coral_cover = centroids_data$high_coral_cover,
    geometry = st_sfc(grid_polys, crs = 4326)
  )
  
  cat("Grid polygons created\n\n")
  
  ################################## Load MATLAB connectivity results ##################################
  
  cat("Loading connectivity matrices from MATLAB...\n")
  
  # Define path to MATLAB cache file (adjust as needed)
  matlab_cache_file <- "G:/Other computers/My MacBook Pro/Documents/Farmer_Ben_Dissertation/Chapters/Chapter3_Metacommunities/temp/conn_structs_20181201_to_20191231_notraj_2018bf.mat"
  
  if (!file.exists(matlab_cache_file)) {
    stop("MATLAB cache file not found: ", matlab_cache_file, 
         "\nPlease run the MATLAB script with RECALCULATE_CONNECTIVITY = true first.")
  }
  
  # Read HDF5 structure
  h5_contents <- h5ls(matlab_cache_file)
  
  # Read original indices - these map from filtered matrix rows back to full reef list
  original_indices <- as.vector(h5read(matlab_cache_file, "/original_indices"))
  
  cat("Connectivity matrices use:", length(original_indices), "reefs\n")
  cat("Original indices range:", range(original_indices), "\n")
  cat("Full reef dataset has:", nrow(centroids_data), "reefs\n\n")
  
  ################################## Identify non-source reefs ##################################
  
  cat("Identifying reefs that never produced trajectories...\n")
  
  # Track which reefs (in the filtered connectivity matrix space) ever had outgoing connections
  n_filtered_reefs <- length(original_indices)
  filtered_sources_ever <- rep(FALSE, n_filtered_reefs)
  
  # Get list of all sparse matrix groups in /#refs#
  refs_contents <- h5_contents[h5_contents$group == "/#refs#" & h5_contents$otype == "H5I_GROUP", ]
  sparse_matrix_names <- refs_contents$name
  
  cat("Found", length(sparse_matrix_names), "sparse matrices in /#refs#\n")
  
  # Check each sparse matrix
  for (i in seq_along(sparse_matrix_names)) {
    matrix_name <- sparse_matrix_names[i]
    matrix_path <- paste0("/#refs#/", matrix_name)
    
    tryCatch({
      # Read sparse matrix row indices
      ir <- h5read(matlab_cache_file, paste0(matrix_path, "/ir"))
      
      # ir contains 0-based row indices into the FILTERED matrix
      if (length(ir) > 0) {
        row_indices <- ir + 1  # Convert to 1-based R indexing
        # Mark these rows in the filtered space as having trajectories
        filtered_sources_ever[unique(row_indices)] <- TRUE
      }
      
    }, error = function(e) {
      # Skip if not a sparse matrix
    })
    
    if (i %% 10 == 0) cat("  Processed", i, "of", length(sparse_matrix_names), "matrices\n")
  }
  
  cat("  Processed all", length(sparse_matrix_names), "matrices\n")
  
  # Map back to full reef dataset
  # Initialize: all reefs assumed to have NO trajectories
  n_reefs <- nrow(centroids_data)
  all_sources_ever <- rep(FALSE, n_reefs)
  
  # Mark reefs that ARE in the filtered set and DO have trajectories
  reefs_in_filtered_set <- original_indices
  all_sources_ever[reefs_in_filtered_set] <- filtered_sources_ever
  
  cat("\nDiagnostic info:\n")
  cat("  Total reefs in full dataset:", n_reefs, "\n")
  cat("  Reefs in connectivity matrices:", n_filtered_reefs, "\n")
  cat("  Of those, reefs with trajectories:", sum(filtered_sources_ever), "\n")
  cat("  Overall, reefs marked as sources:", sum(all_sources_ever), "\n")
  
  # Identify inactive sources
  inactive_sources <- which(!all_sources_ever)
  active_sources <- which(all_sources_ever)
  
  cat("\n=== TRAJECTORY ANALYSIS RESULTS ===\n")
  cat("Reefs that NEVER produced trajectories:", length(inactive_sources), 
      sprintf("(%.1f%%)\n", 100 * length(inactive_sources) / n_reefs))
  cat("Reefs with trajectories:", length(active_sources), 
      sprintf("(%.1f%%)\n", 100 * length(active_sources) / n_reefs))
  
  if (length(inactive_sources) > 0) {
    cat("\nSample inactive reef IDs (first 10):", 
        paste(centroids_data$unique_ID[inactive_sources[1:min(10, length(inactive_sources))]], 
              collapse = ", "), "\n")
  }
  cat("====================================\n\n")
  
  ################################## Filter data ##################################
  
  cat("\nFiltering data to active sources only...\n")
  
  centroids_filtered <- centroids_data[active_sources, ]
  grid_filtered_sf <- grid_sf[active_sources, ]
  grid_inactive_sf <- grid_sf[inactive_sources, ]
  
  cat("Filtered reef sites:", nrow(centroids_filtered), "\n")
  cat("Reduction:", round((1 - nrow(centroids_filtered)/nrow(centroids_data)) * 100, 1), "%\n\n")
  
  ################################## Convert to spatial objects ##################################
  
  cat("Creating spatial objects for mapping...\n")
  
  # Add trajectory status to full grid
  grid_sf$has_trajectories <- all_sources_ever
  
  cat("Spatial objects created\n\n")
  
  ################################## Create leaflet maps ##################################
  
  cat("Creating interactive leaflet maps...\n\n")
  
  # Create color palettes with viridis
  pal_total <- colorNumeric(palette = "viridis", domain = grid_filtered_sf$mean_coral_cover * 100)
  pal_low <- colorNumeric(palette = "plasma", domain = grid_filtered_sf$low_coral_cover * 100)
  pal_moderate <- colorNumeric(palette = "inferno", domain = grid_filtered_sf$moderate_coral_cover * 100)
  pal_high <- colorNumeric(palette = "magma", domain = grid_filtered_sf$high_coral_cover * 100)
  
  # Map 1: Total coral cover
  cat("Map 1: Total coral cover\n")
  map1 <- leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addTiles(group = "OpenStreetMap") %>%
    addPolygons(data = grid_filtered_sf,
                fillColor = ~pal_total(mean_coral_cover * 100),
                fillOpacity = 0.6,
                color = "black",
                weight = 1,
                popup = ~paste0("<b>ID:</b> ", unique_ID, 
                                "<br><b>Total Cover:</b> ", round(mean_coral_cover * 100, 2), "%")) %>%
    addLegend("bottomright", pal = pal_total, 
              values = grid_filtered_sf$mean_coral_cover * 100,
              title = "Total Coral<br>Cover (%)",
              opacity = 1) %>%
    addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
                     options = layersControlOptions(collapsed = FALSE))
  print(map1)
  
  # Map 2: Low susceptibility cover
  cat("Map 2: Low susceptibility cover\n")
  map2 <- leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addTiles(group = "OpenStreetMap") %>%
    addPolygons(data = grid_filtered_sf,
                fillColor = ~pal_low(low_coral_cover * 100),
                fillOpacity = 0.6,
                color = "black",
                weight = 1,
                popup = ~paste0("<b>ID:</b> ", unique_ID, 
                                "<br><b>Low Susceptibility:</b> ", round(low_coral_cover * 100, 2), "%")) %>%
    addLegend("bottomright", pal = pal_low, 
              values = grid_filtered_sf$low_coral_cover * 100,
              title = "Low Susceptibility<br>Cover (%)",
              opacity = 1) %>%
    addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
                     options = layersControlOptions(collapsed = FALSE))
  print(map2)
  
  # Map 3: Moderate susceptibility cover
  cat("Map 3: Moderate susceptibility cover\n")
  map3 <- leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addTiles(group = "OpenStreetMap") %>%
    addPolygons(data = grid_filtered_sf,
                fillColor = ~pal_moderate(moderate_coral_cover * 100),
                fillOpacity = 0.6,
                color = "black",
                weight = 1,
                popup = ~paste0("<b>ID:</b> ", unique_ID, 
                                "<br><b>Moderate Susceptibility:</b> ", round(moderate_coral_cover * 100, 2), "%")) %>%
    addLegend("bottomright", pal = pal_moderate, 
              values = grid_filtered_sf$moderate_coral_cover * 100,
              title = "Moderate Susceptibility<br>Cover (%)",
              opacity = 1) %>%
    addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
                     options = layersControlOptions(collapsed = FALSE))
  print(map3)
  
  # Map 4: High susceptibility cover
  cat("Map 4: High susceptibility cover\n")
  map4 <- leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addTiles(group = "OpenStreetMap") %>%
    addPolygons(data = grid_filtered_sf,
                fillColor = ~pal_high(high_coral_cover * 100),
                fillOpacity = 0.6,
                color = "black",
                weight = 1,
                popup = ~paste0("<b>ID:</b> ", unique_ID, 
                                "<br><b>High Susceptibility:</b> ", round(high_coral_cover * 100, 2), "%")) %>%
    addLegend("bottomright", pal = pal_high, 
              values = grid_filtered_sf$high_coral_cover * 100,
              title = "High Susceptibility<br>Cover (%)",
              opacity = 1) %>%
    addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
                     options = layersControlOptions(collapsed = FALSE))
  print(map4)
  
  # Map 5: Comparison - Active vs. Inactive sources
  cat("Map 5: Active vs. Inactive source reefs\n")
  map5 <- leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
    addTiles(group = "OpenStreetMap") %>%
    addPolygons(data = grid_inactive_sf,
                fillColor = "red",
                fillOpacity = 0.6,
                color = "darkred",
                weight = 1,
                popup = ~paste0("<b>ID:</b> ", unique_ID, 
                                "<br><b>Status:</b> No trajectories")) %>%
    addPolygons(data = grid_filtered_sf,
                fillColor = ~pal_total(mean_coral_cover * 100),
                fillOpacity = 0.6,
                color = "black",
                weight = 1,
                popup = ~paste0("<b>ID:</b> ", unique_ID, 
                                "<br><b>Total Cover:</b> ", round(mean_coral_cover * 100, 2), "%",
                                "<br><b>Status:</b> Active source")) %>%
    addLegend("bottomright", 
              colors = c("red", "blue"),
              labels = c(sprintf("No trajectories (n=%d)", length(inactive_sources)),
                         sprintf("Active sources (n=%d)", nrow(grid_filtered_sf))),
              title = "Reef Status",
              opacity = 1) %>%
    addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
                     options = layersControlOptions(collapsed = FALSE))
  print(map5)
  
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
  
  cat("\n=== COMPLETE ===\n")