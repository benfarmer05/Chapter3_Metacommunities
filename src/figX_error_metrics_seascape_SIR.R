  
  # .rs.restartR(clean = TRUE)
  
  rm(list=ls())
  
  library(here)
  library(rhdf5)
  library(terra)
  library(tidyverse)
  library(sf)
  library(leaflet)
  library(patchwork)
  
  
  # NOTE - 20 Nov 2025
  #   something to fix still:
  #     - some wonky nearest distances going up to northern reefs! not right

  ################################## Settings ##################################
  
  use_cached_distances <- TRUE
  disease_thresh <- 0.001
  obs_extent_lonlat <- c(-65.58409, -64.50770, 18.05, 18.64853)
  
  # Scoring parameters
  model_res = 650
  lambda_P <- model_res*2  # 1300m - presences should be close
  lambda_A <- model_res*6  # 1950m - absences should be far
  w1 <- 0.7
  w2 <- 0.3
  
  # Binary thresholds for Score 2
  binary_threshold_P <- lambda_P  # Presence within this = match
  binary_threshold_A <- lambda_A  # Absence beyond this = match
  
  # Auto-detect scenario folders
  seascape_dir <- here("output", "seascape_SIR")
  all_folders <- list.dirs(seascape_dir, full.names = FALSE, recursive = FALSE)
  scenarios <- all_folders[grepl("SCENARIO[0-9]", all_folders, ignore.case = TRUE)]
  scenarios <- sort(scenarios)
  
  cat("Detected", length(scenarios), "scenario folders:\n")
  for (s in scenarios) cat("  -", s, "\n")
  cat("\n")
  
  ################################## Load Data ##################################
  
  cat("Loading observations...\n")
  obs <- read.csv(here("output", "combined_coral_data.csv")) %>%
    mutate(date = as.POSIXct(date, tz = "UTC"),
           month = floor_date(as.Date(date), "month")) %>%
    filter(year(date) %in% c(2018, 2019), presence %in% c("P", "S", "A"))
  
  cat("Loaded", nrow(obs), "observations\n")
  cat("  Presences:", sum(obs$presence %in% c("P", "S")), "\n")
  cat("  Absences:", sum(obs$presence == "A"), "\n\n")
  
  # Load bathymetry
  spatial_metadata <- readRDS(here("data", "spatial_metadata.rds"))
  bathy_full <- readRDS(here("data", "bathy_final.rds"))
  terra::crs(bathy_full) <- spatial_metadata$bathy_final$crs
  
  obs_extent_vect <- vect(ext(obs_extent_lonlat), crs = "EPSG:4326")
  obs_extent_proj <- project(obs_extent_vect, crs(bathy_full))
  bathy_final <- crop(bathy_full, obs_extent_proj)
  
  if (!dir.exists(here("temp"))) dir.create(here("temp"))
  
  ################################## Create Spatial Objects ##################################
  
  obs_sf <- st_as_sf(obs, coords = c("lon", "lat"), crs = 4326) %>%
    st_transform(crs(bathy_final))
  
  ################################## Filter Observations ##################################
  
  cat("\nFiltering observations to study extent...\n")
  obs_coords_lonlat <- st_coordinates(st_transform(obs_sf, 4326))
  within_extent <- (obs_coords_lonlat[, 1] >= obs_extent_lonlat[1] & 
                      obs_coords_lonlat[, 1] <= obs_extent_lonlat[2] & 
                      obs_coords_lonlat[, 2] >= obs_extent_lonlat[3] & 
                      obs_coords_lonlat[, 2] <= obs_extent_lonlat[4])
  
  n_outside <- sum(!within_extent)
  if (n_outside > 0) cat("  Dropping", n_outside, "observations outside study extent\n")
  
  obs_sf <- obs_sf[within_extent, ]
  obs <- obs[within_extent, ]
  
  cat("  Retained", nrow(obs_sf), "observations\n\n")
  
  ################################## Process All Scenarios ##################################
  
  all_monthly_summaries <- list()
  all_map_data <- list()
  all_diseased_sites <- list()
  
  for (scenario in scenarios) {
    cat("\n========================================\n")
    cat("Processing scenario:", scenario, "\n")
    cat("========================================\n\n")
    
    matlab_file <- here("output", "seascape_SIR", scenario, "seascape_SIR_workspace.mat")
    
    if (!file.exists(matlab_file)) {
      cat("WARNING: File not found, skipping:", matlab_file, "\n")
      next
    }
    
    # Load scenario data
    locations <- h5read(matlab_file, "/outputs/sites/locations")
    if (ncol(locations) > nrow(locations)) locations <- t(locations)
    unique_IDs <- as.vector(h5read(matlab_file, "/outputs/sites/unique_IDs"))
    R_total <- h5read(matlab_file, "/outputs/totals/R_total")
    
    ref_date <- as.Date("2018-12-01")
    tspan <- as.vector(h5read(matlab_file, "/outputs/metadata/tspan"))
    sim_dates <- ref_date + (tspan[1]:tspan[2]) - 1
    
    sites_sf <- st_as_sf(data.frame(unique_ID = unique_IDs, lon = locations[,1], lat = locations[,2]),
                         coords = c("lon", "lat"), crs = 4326) %>%
      st_transform(crs(bathy_final))
    
    unique_months <- sort(unique(obs$month))
    
    monthly_results <- lapply(unique_months, function(month_date) {
      
      cache_file <- here("temp", paste0("distances_", gsub("[^A-Za-z0-9]", "_", scenario), 
                                        "_", format(month_date, "%Y%m"), ".rds"))
      
      obs_month <- obs_sf %>% filter(month == month_date)
      
      month_end <- ceiling_date(month_date, "month") - days(1)
      day_idx <- which.min(abs(sim_dates - month_end))
      diseased_status <- R_total[day_idx, ] >= disease_thresh
      diseased_sites <- sites_sf[diseased_status, ]
      
      if (nrow(diseased_sites) == 0) {
        return(list(
          summary = data.frame(month = month_date, mean_distance_m = NA, 
                               n_obs = nrow(obs_month), n_diseased_sites = 0,
                               n_presence = sum(obs_month$presence %in% c("P", "S")),
                               n_absence = sum(obs_month$presence == "A")),
          map_data = NULL,
          diseased_sites = NULL
        ))
      }
      
      cat("Processing", nrow(obs_month), "observations for", format(month_date, "%b %Y"), 
          "(", sum(obs_month$presence %in% c("P", "S")), "P,", sum(obs_month$presence == "A"), "A)...\n")
      
      if (use_cached_distances && file.exists(cache_file)) {
        cat("  Loading cached distances...\n")
        cached <- readRDS(cache_file)
        nearest_data <- cached$nearest_data
        
      } else {
        source_raster <- bathy_final
        source_raster[!is.na(source_raster)] <- 1
        diseased_cells <- cellFromXY(source_raster, st_coordinates(diseased_sites))
        source_raster[diseased_cells] <- 0
        dist_raster <- gridDist(source_raster, target = 0)
        
        obs_coords <- st_coordinates(obs_month)
        diseased_coords <- st_coordinates(diseased_sites)
        
        cat("    Finding water-nearest diseased site for each observation...\n")
        start_time <- Sys.time()
        
        nearest_data <- map_dfr(1:nrow(obs_month), function(i) {
          obs_source <- bathy_final
          obs_source[!is.na(obs_source)] <- 1
          obs_cell <- cellFromXY(obs_source, obs_coords[i, , drop = FALSE])
          
          if (is.na(obs_source[obs_cell])) {
            water_cells <- which(!is.na(values(obs_source)))
            water_coords <- xyFromCell(obs_source, water_cells)
            dists_to_water <- sqrt(rowSums((t(water_coords) - obs_coords[i,])^2))
            obs_cell <- water_cells[which.min(dists_to_water)]
          }
          
          obs_source[obs_cell] <- 0
          obs_dist_raster <- gridDist(obs_source, target = 0)
          dists_to_diseased <- terra::extract(obs_dist_raster, diseased_sites)[,2]
          
          if (all(is.na(dists_to_diseased))) return(NULL)
          
          nearest_idx <- which.min(dists_to_diseased)
          nearest_dist <- dists_to_diseased[nearest_idx]
          
          data.frame(
            obs_idx = i,
            obs_lon = obs_coords[i, 1],
            obs_lat = obs_coords[i, 2],
            nearest_lon = diseased_coords[nearest_idx, 1],
            nearest_lat = diseased_coords[nearest_idx, 2],
            distance_m = nearest_dist,
            presence = as.character(obs_month$presence[i])
          )
        })
        
        n_disconnected <- nrow(obs_month) - nrow(nearest_data)
        if (n_disconnected > 0) {
          cat("    Skipped", n_disconnected, "observations in isolated water bodies\n")
        }
        
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        cat("    Completed in", round(elapsed, 1), "seconds\n")
        
        # Outlier correction
        n_corrected <- 0
        if (nrow(nearest_data) > 1) {
          for (i in 1:nrow(nearest_data)) {
            dists_between_obs <- sqrt((nearest_data$obs_lon - nearest_data$obs_lon[i])^2 + 
                                        (nearest_data$obs_lat - nearest_data$obs_lat[i])^2)
            nearby_idx <- which(dists_between_obs < 100 & dists_between_obs > 0)
            
            if (length(nearby_idx) > 0) {
              nearby_dists <- nearest_data$distance_m[nearby_idx]
              my_dist <- nearest_data$distance_m[i]
              median_nearby <- median(nearby_dists, na.rm = TRUE)
              
              if (!is.na(median_nearby) && median_nearby > 0 && my_dist > 2 * median_nearby) {
                median_obs_idx <- nearby_idx[which.min(abs(nearby_dists - median_nearby))]
                nearest_data$distance_m[i] <- nearest_data$distance_m[median_obs_idx]
                nearest_data$nearest_lon[i] <- nearest_data$nearest_lon[median_obs_idx]
                nearest_data$nearest_lat[i] <- nearest_data$nearest_lat[median_obs_idx]
                n_corrected <- n_corrected + 1
              }
            }
          }
          if (n_corrected > 0) {
            cat("    Corrected", n_corrected, "suspicious distances\n")
          }
        }
        
        cat("    Saving to cache:", basename(cache_file), "\n")
        saveRDS(list(nearest_data = nearest_data), cache_file)
      }
      
      obs_wgs84 <- st_transform(obs_month, 4326)
      diseased_wgs84 <- st_transform(diseased_sites, 4326)
      nearest_wgs84 <- st_as_sf(nearest_data[, c("nearest_lon", "nearest_lat")], 
                                coords = c("nearest_lon", "nearest_lat"), 
                                crs = crs(bathy_final)) %>%
        st_transform(4326)
      
      map_data <- data.frame(
        month = month_date,
        month_label = format(month_date, "%b %Y"),
        obs_lon = st_coordinates(obs_wgs84)[,1],
        obs_lat = st_coordinates(obs_wgs84)[,2],
        nearest_lon = st_coordinates(nearest_wgs84)[,1],
        nearest_lat = st_coordinates(nearest_wgs84)[,2],
        distance_m = nearest_data$distance_m,
        presence = nearest_data$presence
      )
      
      diseased_coords_wgs84 <- st_coordinates(diseased_wgs84)
      
      list(
        summary = data.frame(
          month = month_date,
          mean_distance_m = mean(nearest_data$distance_m, na.rm = TRUE),
          n_obs = nrow(obs_month),
          n_diseased_sites = nrow(diseased_sites),
          n_presence = sum(obs_month$presence %in% c("P", "S")),
          n_absence = sum(obs_month$presence == "A"),
          mean_dist_presence = mean(nearest_data$distance_m[nearest_data$presence %in% c("P", "S")], na.rm = TRUE),
          mean_dist_absence = mean(nearest_data$distance_m[nearest_data$presence == "A"], na.rm = TRUE),
          # Score 1: Continuous distance-based
          score_P = mean(exp(-nearest_data$distance_m[nearest_data$presence %in% c("P", "S")] / lambda_P), na.rm = TRUE),
          score_A = mean(1 - 1 / (1 + nearest_data$distance_m[nearest_data$presence == "A"] / lambda_A), na.rm = TRUE),
          score_combined = w1 * mean(exp(-nearest_data$distance_m[nearest_data$presence %in% c("P", "S")] / lambda_P), na.rm = TRUE) +
            w2 * mean(1 - 1 / (1 + nearest_data$distance_m[nearest_data$presence == "A"] / lambda_A), na.rm = TRUE),
          # Score 2: Binary threshold-based
          n_presence_match = sum(nearest_data$distance_m[nearest_data$presence %in% c("P", "S")] <= binary_threshold_P, na.rm = TRUE),
          n_absence_match = sum(nearest_data$distance_m[nearest_data$presence == "A"] >= binary_threshold_A, na.rm = TRUE),
          alpha = sum(nearest_data$presence %in% c("P", "S")) / sum(nearest_data$presence == "A"),
          score2_raw = sum(nearest_data$distance_m[nearest_data$presence %in% c("P", "S")] <= binary_threshold_P, na.rm = TRUE) +
            (sum(nearest_data$presence %in% c("P", "S")) / sum(nearest_data$presence == "A")) * 
            sum(nearest_data$distance_m[nearest_data$presence == "A"] >= binary_threshold_A, na.rm = TRUE),
          score2_max = sum(nearest_data$presence %in% c("P", "S")) + 
            (sum(nearest_data$presence %in% c("P", "S")) / sum(nearest_data$presence == "A")) * 
            sum(nearest_data$presence == "A"),
          score2_norm = (sum(nearest_data$distance_m[nearest_data$presence %in% c("P", "S")] <= binary_threshold_P, na.rm = TRUE) +
                           (sum(nearest_data$presence %in% c("P", "S")) / sum(nearest_data$presence == "A")) * 
                           sum(nearest_data$distance_m[nearest_data$presence == "A"] >= binary_threshold_A, na.rm = TRUE)) /
            (sum(nearest_data$presence %in% c("P", "S")) + 
               (sum(nearest_data$presence %in% c("P", "S")) / sum(nearest_data$presence == "A")) * 
               sum(nearest_data$presence == "A"))
        ),
        map_data = map_data,
        diseased_sites = data.frame(
          month_label = format(month_date, "%b %Y"),
          lon = diseased_coords_wgs84[,1],
          lat = diseased_coords_wgs84[,2]
        )
      )
    })
    
    monthly_summary <- bind_rows(lapply(monthly_results, function(x) x$summary))
    monthly_summary$scenario <- scenario
    
    all_monthly_summaries[[scenario]] <- monthly_summary
    all_map_data[[scenario]] <- bind_rows(lapply(monthly_results, function(x) x$map_data))
    all_diseased_sites[[scenario]] <- bind_rows(lapply(monthly_results, function(x) x$diseased_sites))
    
    cat("\n=== MONTHLY SUMMARY:", scenario, "===\n")
    print(monthly_summary)
  }
  
  combined_summary <- bind_rows(all_monthly_summaries)
  
  ################################## Plot Results ##################################
  
  # Plot 1: Mean distances by presence/absence for all scenarios
  distance_plot_split <- combined_summary %>%
    select(scenario, month, mean_dist_presence, mean_dist_absence, n_presence, n_absence) %>%
    pivot_longer(cols = c(mean_dist_presence, mean_dist_absence), 
                 names_to = "type", values_to = "distance") %>%
    mutate(type = ifelse(type == "mean_dist_presence", "Presence", "Absence")) %>%
    ggplot(aes(x = month, y = distance, color = type, linetype = scenario, group = interaction(type, scenario))) +
    geom_line(linewidth = 1) +
    geom_point(alpha = 0.5, size = 2) +
    scale_color_manual(values = c("Presence" = "coral", "Absence" = "steelblue"),
                       name = "Observation Type") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Scenario") +
    scale_y_continuous(labels = scales::comma) +
    labs(title = "Mean Water Distance: Observations to Nearest Diseased Model Site",
         subtitle = "Disease threshold: 0.1% cover lost",
         x = "Month", y = "Mean Distance (m)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"),
          legend.position = "right")
  
  print(distance_plot_split)
  
  # Score plots
  score_P_plot <- ggplot(combined_summary, aes(x = month, y = score_P, color = scenario, linetype = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(alpha = 0.5, size = 2) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Presence Score (Score_P)",
         subtitle = paste0("Exponential decay: exp(-distance/", lambda_P, ")"),
         x = "Month", y = "Score (0-1, higher = better)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
  
  score_A_plot <- ggplot(combined_summary, aes(x = month, y = score_A, color = scenario, linetype = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(alpha = 0.5, size = 2) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Absence Score (Score_A)",
         subtitle = paste0("Flipped inverse: 1 - 1/(1+distance/", lambda_A, ")"),
         x = "Month", y = "Score (0-1, higher = better)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
  
  combined_scores <- score_P_plot / score_A_plot + 
    plot_annotation(title = "Score 1: Continuous Distance-Based Performance",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  
  print(combined_scores)
  
  # Combined score 1
  combined_score_plot <- ggplot(combined_summary, aes(x = month, y = score_combined, color = scenario, linetype = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(alpha = 0.5, size = 2) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Combined Score 1",
         subtitle = paste0("Weighted: ", w1, " × Score_P + ", w2, " × Score_A"),
         x = "Month", y = "Score (0-1, higher = better)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
  
  print(combined_score_plot)
  
  # Score 2: Binary threshold-based
  score2_plot <- ggplot(combined_summary, aes(x = month, y = score2_norm, color = scenario, linetype = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(alpha = 0.5, size = 2) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Score 2: Binary Threshold Performance",
         subtitle = paste0("Presence ≤", binary_threshold_P, "m, Absence ≥", binary_threshold_A, "m"),
         x = "Month", y = "Normalized Score (0-1, higher = better)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"))
  
  print(score2_plot)
  
  ################################## Leaflet Maps (One per Scenario) ##################################
  
  for (scenario in scenarios) {
    cat("\nCreating leaflet map for:", scenario, "\n")
    
    monthly_map_data <- all_map_data[[scenario]] %>% filter(!is.na(distance_m))
    diseased_sites <- all_diseased_sites[[scenario]] %>% filter(!is.na(lon))
    
    unique_labels <- format(sort(unique(monthly_map_data$month)), "%b %Y")
    
    leaflet_map <- leaflet() %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
      addTiles(group = "OpenStreetMap")
    
    for (month_label in unique_labels) {
      month_data <- monthly_map_data %>% filter(month_label == !!month_label)
      month_diseased <- diseased_sites %>% filter(month_label == !!month_label)
      
      month_presence <- month_data %>% filter(presence %in% c("P", "S"))
      month_absence <- month_data %>% filter(presence == "A")
      
      leaflet_map <- leaflet_map %>%
        addCircleMarkers(
          data = month_diseased,
          lng = ~lon, lat = ~lat,
          radius = 3,
          color = "gray",
          fillColor = "gray",
          fillOpacity = 0.4,
          stroke = TRUE,
          weight = 1,
          group = month_label,
          popup = ~paste0("<b>Diseased Model Site</b><br><b>Month:</b> ", month_label)
        )
      
      if (nrow(month_presence) > 0) {
        leaflet_map <- leaflet_map %>%
          addCircleMarkers(
            data = month_presence,
            lng = ~obs_lon, lat = ~obs_lat,
            radius = 7,
            color = "darkorange",
            fillColor = "coral",
            fillOpacity = 0.9,
            stroke = TRUE,
            weight = 2,
            group = month_label,
            popup = ~paste0("<b>Presence Observation</b><br>",
                            "<b>Month:</b> ", month_label, "<br>",
                            "<b>Water distance:</b> ", round(distance_m), " m")
          )
        
        for (i in 1:nrow(month_presence)) {
          leaflet_map <- leaflet_map %>%
            addPolylines(
              lng = c(month_presence$obs_lon[i], month_presence$nearest_lon[i]),
              lat = c(month_presence$obs_lat[i], month_presence$nearest_lat[i]),
              color = "orange",
              weight = 2,
              opacity = 0.7,
              group = month_label
            )
        }
      }
      
      if (nrow(month_absence) > 0) {
        leaflet_map <- leaflet_map %>%
          addCircleMarkers(
            data = month_absence,
            lng = ~obs_lon, lat = ~obs_lat,
            radius = 7,
            color = "darkblue",
            fillColor = "lightblue",
            fillOpacity = 0.9,
            stroke = TRUE,
            weight = 2,
            group = month_label,
            popup = ~paste0("<b>Absence Observation</b><br>",
                            "<b>Month:</b> ", month_label, "<br>",
                            "<b>Water distance:</b> ", round(distance_m), " m")
          )
        
        for (i in 1:nrow(month_absence)) {
          leaflet_map <- leaflet_map %>%
            addPolylines(
              lng = c(month_absence$obs_lon[i], month_absence$nearest_lon[i]),
              lat = c(month_absence$obs_lat[i], month_absence$nearest_lat[i]),
              color = "blue",
              weight = 2,
              opacity = 0.7,
              group = month_label
            )
        }
      }
    }
    
    leaflet_map <- leaflet_map %>%
      addLegend("bottomright",
                colors = c("coral", "lightblue", "gray"),
                labels = c("Presence Observation", "Absence Observation", 
                           "Diseased Model Sites"),
                title = paste("Scenario:", scenario),
                opacity = 1) %>%
      addLayersControl(
        baseGroups = c("Satellite", "OpenStreetMap"),
        overlayGroups = unique_labels,
        options = layersControlOptions(collapsed = FALSE)
      ) %>%
      hideGroup(unique_labels[-length(unique_labels)])
    
    print(leaflet_map)
  }
  
  cat("\n✓ Script complete!\n")