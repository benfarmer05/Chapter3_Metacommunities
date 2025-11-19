# .rs.restartR(clean = TRUE)

rm(list=ls())

library(here)
library(rhdf5)
library(terra)
library(tidyverse)
library(sf)
library(leaflet)

################################## Settings ##################################

use_cached_distances <- TRUE  # Set TRUE to load cached distances, FALSE to recalculate
disease_thresh <- 0.001

# Spatial extent filter for observations (lon: xmin, xmax; lat: ymin, ymax)
obs_extent_lonlat <- c(-65.58409, -64.50770, 18.05, 18.64853)

# Scoring parameters
model_res = 650 #in meters
lambda_P <- model_res*2  # Decay rate for presences (meters)
lambda_A <- model_res*3  # Decay rate for absences (meters)
w1 <- 0.5  # Weight for presence score
w2 <- 0.5  # Weight for absence score

################################## Load Data ##################################

cat("Loading MATLAB outputs...\n")
matlab_file <- here("output", "seascape_SIR", "seascape_SIR_workspace.mat")

# Extract model grid and disease time series
locations <- h5read(matlab_file, "/outputs/sites/locations")
if (ncol(locations) > nrow(locations)) locations <- t(locations)
unique_IDs <- as.vector(h5read(matlab_file, "/outputs/sites/unique_IDs"))
R_total <- h5read(matlab_file, "/outputs/totals/R_total")

# Create date range
ref_date <- as.Date("2018-12-01")
tspan <- as.vector(h5read(matlab_file, "/outputs/metadata/tspan"))
sim_dates <- ref_date + (tspan[1]:tspan[2]) - 1

# Load observations (include both presences AND absences)
obs <- read.csv(here("output", "combined_coral_data.csv")) %>%
  mutate(date = as.POSIXct(date, tz = "UTC"),
         month = floor_date(as.Date(date), "month")) %>%
  filter(year(date) %in% c(2018, 2019), presence %in% c("P", "S", "A"))

cat("Loaded", ncol(R_total), "sites,", nrow(R_total), "days,", nrow(obs), "observations\n")
cat("  Presences:", sum(obs$presence %in% c("P", "S")), "\n")
cat("  Absences:", sum(obs$presence == "A"), "\n\n")

# Load and prepare bathymetry
spatial_metadata <- readRDS(here("data", "spatial_metadata.rds"))
bathy_full <- readRDS(here("data", "bathy_final.rds"))
terra::crs(bathy_full) <- spatial_metadata$bathy_final$crs

# Crop bathymetry to study extent
obs_extent_vect <- vect(ext(obs_extent_lonlat), crs = "EPSG:4326")
obs_extent_proj <- project(obs_extent_vect, crs(bathy_full))
bathy_final <- crop(bathy_full, obs_extent_proj)

cat("Bathymetry cropped to study extent:", 
    round(ncell(bathy_final) / 1e6, 2), "M cells\n")

# Create temp directory for caching if needed
if (!dir.exists(here("temp"))) {
  dir.create(here("temp"))
}

################################## Create Spatial Objects ##################################

sites_sf <- st_as_sf(data.frame(unique_ID = unique_IDs, lon = locations[,1], lat = locations[,2]),
                     coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs(bathy_final))

obs_sf <- st_as_sf(obs, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs(bathy_final))

################################## Filter Out Invalid Observations ##################################

cat("\nFiltering observations to study extent...\n")

# Get observation coordinates in lon/lat
obs_coords_lonlat <- st_coordinates(st_transform(obs_sf, 4326))

# Check which observations are within specified extent
within_extent <- (obs_coords_lonlat[, 1] >= obs_extent_lonlat[1] & 
                    obs_coords_lonlat[, 1] <= obs_extent_lonlat[2] & 
                    obs_coords_lonlat[, 2] >= obs_extent_lonlat[3] & 
                    obs_coords_lonlat[, 2] <= obs_extent_lonlat[4])

n_outside <- sum(!within_extent)
if (n_outside > 0) {
  cat("  Dropping", n_outside, "observations outside study extent\n")
}

# Keep only observations within extent
obs_sf <- obs_sf[within_extent, ]
obs <- obs[within_extent, ]

cat("  Retained", nrow(obs_sf), "observations within study area\n")
cat("    Presences:", sum(obs$presence %in% c("P", "S")), "\n")
cat("    Absences:", sum(obs$presence == "A"), "\n\n")

################################## Calculate Monthly Distances ##################################

cat("Calculating monthly distances...\n")

disease_thresh <- 0.001
unique_months <- sort(unique(obs$month))

monthly_results <- lapply(unique_months, function(month_date) {
  
  cache_file <- here("temp", paste0("distances_", format(month_date, "%Y%m"), ".rds"))
  
  # Get observations for this month
  obs_month <- obs_sf %>% filter(month == month_date)
  
  # Find diseased sites (use last day of month)
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
  
  # Check cache
  if (use_cached_distances && file.exists(cache_file)) {
    cat("  Loading cached distances from", basename(cache_file), "...\n")
    cached <- readRDS(cache_file)
    distances_m <- cached$distances_m
    nearest_data <- cached$nearest_data
    
  } else {
    # Create distance raster from all diseased sites
    source_raster <- bathy_final
    source_raster[!is.na(source_raster)] <- 1
    diseased_cells <- cellFromXY(source_raster, st_coordinates(diseased_sites))
    source_raster[diseased_cells] <- 0
    
    dist_raster <- gridDist(source_raster, target = 0)
    
    # Extract distances at observation locations (to nearest diseased site by water)
    distances_m <- terra::extract(dist_raster, obs_month)[,2]
    
    # For each observation, find which diseased site is water-distance nearest
    obs_coords <- st_coordinates(obs_month)
    diseased_coords <- st_coordinates(diseased_sites)
    
    cat("    Finding water-nearest diseased site for each observation (est.", 
        round(nrow(obs_month) * 0.5), "sec)...\n")
    
    start_time <- Sys.time()
    
    nearest_data <- map_dfr(1:nrow(obs_month), function(i) {
      # Create distance raster starting from this observation
      obs_source <- bathy_final
      obs_source[!is.na(obs_source)] <- 1
      obs_cell <- cellFromXY(obs_source, obs_coords[i, , drop = FALSE])
      
      # Check if observation is on land (snap to nearest water cell)
      if (is.na(obs_source[obs_cell])) {
        water_cells <- which(!is.na(values(obs_source)))
        water_coords <- xyFromCell(obs_source, water_cells)
        dists_to_water <- sqrt(rowSums((t(water_coords) - obs_coords[i,])^2))
        nearest_water_cell <- water_cells[which.min(dists_to_water)]
        obs_cell <- nearest_water_cell
      }
      
      obs_source[obs_cell] <- 0
      obs_dist_raster <- gridDist(obs_source, target = 0)
      
      # Extract distances at all diseased sites
      dists_to_diseased <- terra::extract(obs_dist_raster, diseased_sites)[,2]
      
      # Skip if disconnected from all diseased sites
      if (all(is.na(dists_to_diseased))) {
        return(NULL)
      }
      
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
    
    # Check for suspicious distance outliers among nearby observations
    n_corrected <- 0
    if (nrow(nearest_data) > 1) {
      cat("    Checking for distance outliers among", nrow(nearest_data), "observations...\n")
      for (i in 1:nrow(nearest_data)) {
        # Find observations within 100m
        dists_between_obs <- sqrt((nearest_data$obs_lon - nearest_data$obs_lon[i])^2 + 
                                    (nearest_data$obs_lat - nearest_data$obs_lat[i])^2)
        nearby_idx <- which(dists_between_obs < 100 & dists_between_obs > 0)
        
        if (length(nearby_idx) > 0) {
          nearby_dists <- nearest_data$distance_m[nearby_idx]
          my_dist <- nearest_data$distance_m[i]
          median_nearby <- median(nearby_dists, na.rm = TRUE)
          
          # Debug for observations with nearby neighbors
          if (length(nearby_idx) >= 3 || my_dist > 10000) {
            cat("      Obs", i, ": dist =", round(my_dist), "m,", 
                length(nearby_idx), "nearby, median =", round(median_nearby), "m\n")
          }
          
          # If this obs has a distance >2x the median of nearby obs, flag it
          if (!is.na(median_nearby) && median_nearby > 0 && my_dist > 2 * median_nearby) {
            # Find which nearby obs has the median distance
            median_obs_idx <- nearby_idx[which.min(abs(nearby_dists - median_nearby))]
            
            cat("      >>> CORRECTING Obs", i, "from", round(my_dist), 
                "m to", round(median_nearby), "m (using nearest site from obs", 
                median_obs_idx, ") <<<\n")
            
            # Copy distance AND nearest site coordinates from the nearby observation
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
    
    # Save to cache
    cat("    Saving to cache:", basename(cache_file), "\n")
    saveRDS(list(distances_m = distances_m, nearest_data = nearest_data), cache_file)
  }
  
  # Convert to WGS84 for leaflet
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
  
  # Store all diseased sites for this month
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
      # Score metrics
      score_P = mean(exp(-nearest_data$distance_m[nearest_data$presence %in% c("P", "S")] / lambda_P), na.rm = TRUE),
      score_A = mean(1 / (1 + nearest_data$distance_m[nearest_data$presence == "A"] / lambda_A), na.rm = TRUE),
      score_combined = w1 * mean(exp(-nearest_data$distance_m[nearest_data$presence %in% c("P", "S")] / lambda_P), na.rm = TRUE) +
        w2 * mean(1 / (1 + nearest_data$distance_m[nearest_data$presence == "A"] / lambda_A), na.rm = TRUE)
    ),
    map_data = map_data,
    diseased_sites = data.frame(
      month_label = format(month_date, "%b %Y"),
      lon = diseased_coords_wgs84[,1],
      lat = diseased_coords_wgs84[,2]
    )
  )
})

# Extract summary and map data
monthly_summary <- bind_rows(lapply(monthly_results, function(x) x$summary))
monthly_map_data <- bind_rows(lapply(monthly_results, function(x) x$map_data))
all_diseased_sites <- bind_rows(lapply(monthly_results, function(x) x$diseased_sites))

cat("\n=== MONTHLY DISTANCE SUMMARY ===\n")
print(monthly_summary)
cat("================================\n\n")

################################## Plot Results ##################################

# Plot 1: Mean distances by presence/absence
distance_plot_split <- monthly_summary %>%
  select(month, mean_dist_presence, mean_dist_absence, n_presence, n_absence) %>%
  pivot_longer(cols = c(mean_dist_presence, mean_dist_absence), 
               names_to = "type", values_to = "distance") %>%
  mutate(type = ifelse(type == "mean_dist_presence", "Presence", "Absence"),
         n = ifelse(type == "Presence", n_presence, n_absence)) %>%
  ggplot(aes(x = month, y = distance, color = type, group = type)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = n), alpha = 0.7) +
  scale_color_manual(values = c("Presence" = "coral", "Absence" = "steelblue"),
                     name = "Observation Type") +
  scale_size_continuous(name = "# Observations", range = c(3, 10)) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Mean Water Distance: Observations to Nearest Diseased Model Site",
       subtitle = "Disease threshold: 0.1% cover lost",
       x = "Month", y = "Mean Distance (m)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"),
        legend.position = "right")

print(distance_plot_split)

# Plot 2: Overall mean (combined)
distance_plot <- ggplot(monthly_summary, aes(x = month, y = mean_distance_m)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(aes(size = n_obs), color = "steelblue", alpha = 0.7) +
  scale_size_continuous(name = "# Observations", range = c(3, 10)) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Mean Water Distance: All Observations to Nearest Diseased Model Site",
       subtitle = "Disease threshold: 0.1% cover lost",
       x = "Month", y = "Mean Distance (m)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

print(distance_plot)

################################## Score Plots ##################################

# Create separate score plots
score_P_plot <- ggplot(monthly_summary, aes(x = month, y = score_P)) +
  geom_line(color = "coral", linewidth = 1) +
  geom_point(aes(size = n_presence), color = "coral", alpha = 0.7) +
  scale_size_continuous(name = "# Presences", range = c(3, 10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Presence Score (Score_P)",
       subtitle = paste0("Exponential decay: exp(-distance/", lambda_P, ")"),
       x = "Month", y = "Score (0-1, higher = better)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

score_A_plot <- ggplot(monthly_summary, aes(x = month, y = score_A)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(aes(size = n_absence), color = "steelblue", alpha = 0.7) +
  scale_size_continuous(name = "# Absences", range = c(3, 10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Absence Score (Score_A)",
       subtitle = paste0("Inverse transform: 1/(1+distance/", lambda_A, ")"),
       x = "Month", y = "Score (0-1, higher = better)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

# Combined plot with both scores
library(patchwork)
combined_scores <- score_P_plot / score_A_plot + 
  plot_annotation(title = "Model Performance Scores Over Time",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

print(combined_scores)

# Overall combined score
combined_score_plot <- ggplot(monthly_summary, aes(x = month, y = score_combined)) +
  geom_line(color = "purple", linewidth = 1) +
  geom_point(aes(size = n_obs), color = "purple", alpha = 0.7) +
  scale_size_continuous(name = "# Observations", range = c(3, 10)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Combined Model Score",
       subtitle = paste0("Weighted: ", w1, " × Score_P + ", w2, " × Score_A"),
       x = "Month", y = "Score (0-1, higher = better)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"))

print(combined_score_plot)

cat("\nMean distance across all months:", 
    round(mean(monthly_summary$mean_distance_m, na.rm = TRUE)), "m\n")

################################## Interactive Leaflet Map ##################################

cat("\nCreating interactive leaflet map...\n")

# Filter out any NULL entries
monthly_map_data <- monthly_map_data %>% filter(!is.na(distance_m))
all_diseased_sites <- all_diseased_sites %>% filter(!is.na(lon))

# Get unique month labels (in chronological order)
unique_labels <- format(sort(unique(monthly_map_data$month)), "%b %Y")

# Create base map
leaflet_map <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addTiles(group = "OpenStreetMap")

# Add layers for each month
for (month_label in unique_labels) {
  month_data <- monthly_map_data %>% filter(month_label == !!month_label)
  month_diseased <- all_diseased_sites %>% filter(month_label == !!month_label)
  
  # Split by presence/absence
  month_presence <- month_data %>% filter(presence %in% c("P", "S"))
  month_absence <- month_data %>% filter(presence == "A")
  
  # Add all diseased sites (small gray circles)
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
  
  # Add PRESENCE observation points (coral/orange)
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
    
    # Add lines for presences (orange)
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
  
  # Add ABSENCE observation points (blue)
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
    
    # Add lines for absences (blue)
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

# Add simple legend and controls (no color scale)
leaflet_map <- leaflet_map %>%
  addLegend("bottomright",
            colors = c("coral", "lightblue", "gray"),
            labels = c("Presence Observation", "Absence Observation", 
                       "Diseased Model Sites"),
            title = "Legend",
            opacity = 1) %>%
  addLayersControl(
    baseGroups = c("Satellite", "OpenStreetMap"),
    overlayGroups = unique_labels,
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  hideGroup(unique_labels[-length(unique_labels)])

print(leaflet_map)

cat("\n✓ Script complete!\n")