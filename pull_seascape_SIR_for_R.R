# .rs.restartR(clean = TRUE)
rm(list=ls())

library(here)
library(rhdf5)
library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(viridis)

################################## Load MATLAB outputs (HDF5 format) ##################################

cat("Loading MATLAB workspace (HDF5 format)...\n")
matlab_file <- here("output", "seascape_SIR", "seascape_SIR_workspace.mat")

if (!file.exists(matlab_file)) {
  stop("MATLAB file not found: ", matlab_file)
}

# Read HDF5 structure
h5_contents <- h5ls(matlab_file)

# Extract locations (MATLAB stores as 2 x n_sites, need to transpose)
locations <- h5read(matlab_file, "/outputs/sites/locations")
if (ncol(locations) > nrow(locations)) {
  locations <- t(locations)  # Now rows = sites, cols = [lon, lat]
}

# Extract unique IDs
unique_IDs <- as.vector(h5read(matlab_file, "/outputs/sites/unique_IDs"))

# Extract N_site
N_site <- as.vector(h5read(matlab_file, "/outputs/sites/N_site"))

# Extract disease states (S, I, R) - already in correct format (days x sites)
S_total <- h5read(matlab_file, "/outputs/totals/S_total")
I_total <- h5read(matlab_file, "/outputs/totals/I_total")
R_total <- h5read(matlab_file, "/outputs/totals/R_total")

num_days <- nrow(R_total)
num_sites <- ncol(R_total)

# Extract dates - reconstruct from what we know
# MATLAB ref = datetime('01-Dec-2018'), so we'll use that
ref_date <- as.Date("2018-12-01")

tspan <- as.vector(h5read(matlab_file, "/outputs/metadata/tspan"))
# MATLAB: ref + days(tspan(1):tspan(2)) creates dates starting from tspan(1)
sim_dates <- ref_date + (tspan[1]:tspan[2]) - 1  # -1 because MATLAB days() is 0-indexed

num_days <- nrow(R_total)
num_sites <- ncol(R_total)

cat("Loaded", num_sites, "sites,", num_days, "days\n")
cat("Date range:", min(sim_dates), "to", max(sim_dates), "\n\n")

################################## Load observations ##################################

cat("Loading observations...\n")
obs <- read.csv(here("output", "combined_coral_data.csv")) %>%
  mutate(date = as.POSIXct(date, tz = "UTC"),
         presence = factor(presence))

# Filter to 2019 and presence only
obs_2019 <- obs %>%
  filter(year(date) == 2019, presence %in% c("P", "S")) %>%
  arrange(date)

cat("Loaded", nrow(obs_2019), "presence observations in 2019\n\n")

################################## Create spatial objects ##################################

# Model sites as spatial points (using terra)
sites_vect <- vect(data.frame(
  unique_ID = unique_IDs,
  lon = locations[, 1],
  lat = locations[, 2],
  N_site = N_site
), geom = c("lon", "lat"), crs = "EPSG:4326")

# Observations as spatial points (using terra)
obs_vect <- vect(obs_2019, geom = c("lon", "lat"), crs = "EPSG:4326")

################################## Calculate disease presence at sites over time ##################################

cat("Calculating disease presence at model sites...\n")

# Define threshold for "diseased" - sites with >0.1% cover lost
disease_thresh <- 0.001

# For each time point, identify which sites are diseased
disease_status <- R_total >= disease_thresh

# Find first day each site became diseased (arrival time)
arrival_day <- apply(disease_status, 2, function(site_vals) {
  first_diseased <- which(site_vals)[1]
  if (is.na(first_diseased)) return(NA)
  return(first_diseased)
})

# Convert to dates
arrival_date <- sim_dates[arrival_day]

# Add to sites data
sites_vect$arrival_day <- arrival_day
sites_vect$arrival_date <- arrival_date
sites_vect$ever_diseased <- !is.na(arrival_day)

cat("Sites ever diseased:", sum(sites_vect$ever_diseased), 
    sprintf("(%.1f%%)\n\n", 100 * mean(sites_vect$ever_diseased)))

################################## Match observations to nearest model sites ##################################

cat("Matching observations to model sites...\n")

# Find nearest model site for each observation using terra
nearest_result <- nearest(obs_vect, sites_vect)

# Extract the indices and distances
nearest_idx <- nearest_result$to_id
distances <- nearest_result$distance

obs_matched <- obs_2019 %>%
  mutate(
    nearest_site = unique_IDs[nearest_idx],
    distance_m = as.numeric(distances),
    obs_date = as.Date(date),
    site_arrival_date = arrival_date[nearest_idx],
    site_ever_diseased = !is.na(site_arrival_date),
    # Model prediction at observation time
    days_since_start = as.numeric(obs_date - min(sim_dates)) + 1,
    # Clamp to valid day range
    days_since_start = pmin(pmax(days_since_start, 1), num_days)
  )

# Extract model prediction at observation time
obs_matched$model_diseased <- FALSE
for (i in 1:nrow(obs_matched)) {
  day_idx <- obs_matched$days_since_start[i]
  site_idx <- nearest_idx[i]
  obs_matched$model_diseased[i] <- disease_status[day_idx, site_idx]
}

cat("Observations matched to sites\n")
cat("Mean distance to nearest site:", round(mean(obs_matched$distance_m)), "m\n\n")

################################## Calculate error metrics ##################################

cat("=== ERROR METRICS ===\n")

# Overall accuracy
accuracy <- mean(obs_matched$model_diseased)
cat("Model detected disease:", sum(obs_matched$model_diseased), "of", nrow(obs_matched), 
    sprintf("observations (%.1f%%)\n", 100 * accuracy))

# Temporal bias
obs_matched_both <- obs_matched %>%
  filter(site_ever_diseased & model_diseased)

if (nrow(obs_matched_both) > 0) {
  temporal_error <- as.numeric(obs_matched_both$site_arrival_date - obs_matched_both$obs_date)
  cat("Mean arrival time error:", round(mean(temporal_error, na.rm = TRUE)), "days\n")
  cat("  (positive = model too slow, negative = model too fast)\n")
  cat("Median arrival time error:", round(median(temporal_error, na.rm = TRUE)), "days\n")
}

# Spatial error - distance to nearest diseased site when model missed
obs_missed <- obs_matched %>% filter(!model_diseased)

if (nrow(obs_missed) > 0) {
  cat("\nObservations where model missed disease:", nrow(obs_missed), "\n")
  cat("Mean distance to nearest site:", round(mean(obs_missed$distance_m)), "m\n")
}

cat("=====================\n\n")

################################## Quick visualization ##################################

cat("Creating comparison map...\n")

# Load land data
land_vect <- readRDS(here("output", "osm_land_vect.rds"))

# Load the original 650m grid
centroids_data <- read.csv(here("data", "centroids_vertices_FINALFORCMS.csv"))

# Final day disease status
final_day_idx <- num_days
final_diseased <- disease_status[final_day_idx, ]

sites_vect$final_diseased <- final_diseased
sites_vect$final_R <- R_total[final_day_idx, ]

# Match grid data to model sites by unique_ID
grid_data <- centroids_data[match(sites_vect$unique_ID, centroids_data$unique_ID), ]
grid_data$final_R <- sites_vect$final_R
grid_data$N_site <- sites_vect$N_site

# Calculate relative cover lost (proportion of initial cover)
grid_data$relative_R <- ifelse(grid_data$N_site > 0, 
                               grid_data$final_R / grid_data$N_site, 
                               NA)

# Create polygons from vertices
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

# Remove NULLs and create sf object
valid_idx <- !sapply(grid_polys, is.null)
grid_sf <- st_sf(
  unique_ID = grid_data$unique_ID[valid_idx],
  final_R = grid_data$final_R[valid_idx],
  relative_R = grid_data$relative_R[valid_idx],
  geometry = st_sfc(grid_polys[valid_idx], crs = 4326)
)

# Get observation coordinates
obs_coords <- crds(obs_vect)

# Calculate extent of grid squares
grid_bbox <- st_bbox(grid_sf)
xlim_grid <- c(grid_bbox["xmin"], grid_bbox["xmax"])
ylim_grid <- c(grid_bbox["ymin"], grid_bbox["ymax"])

# Original extent
xlim_orig <- c(-65.39, -64.65)
ylim_orig <- c(18.06, 18.47)

# Calculate max relative loss rounded up to nearest 0.05
max_relative <- max(grid_sf$relative_R, na.rm = TRUE)
max_relative_rounded <- ceiling(max_relative / 0.05) * 0.05

cat("Max relative loss:", round(max_relative, 3), "-> rounded to", max_relative_rounded, "\n")

# Plot 1: Absolute cover lost (original extent)
plot1 <- ggplot() +
  geom_sf(data = grid_sf, aes(fill = final_R), color = NA) +
  scale_fill_viridis(name = "Cover lost\n(absolute)", option = "plasma", 
                     na.value = "transparent") +
  geom_spatvector(data = land_vect, fill = "gray90", color = "black") +
  geom_point(data = as.data.frame(obs_coords), 
             aes(x = x, y = y), 
             shape = 4, size = 2, color = "red", stroke = 1.5) +
  coord_sf(xlim = xlim_orig, ylim = ylim_orig, expand = FALSE) +
  labs(title = paste("Absolute Cover Lost -", max(sim_dates)),
       subtitle = "Red X = observed disease presence") +
  theme_classic(base_family = "Georgia")

# Plot 2: Relative cover lost (original extent)
plot2 <- ggplot() +
  geom_sf(data = grid_sf, aes(fill = relative_R), color = NA) +
  scale_fill_viridis(name = "Cover lost\n(proportion)", option = "plasma", 
                     na.value = "transparent", limits = c(0, max_relative_rounded)) +
  geom_spatvector(data = land_vect, fill = "gray90", color = "black") +
  geom_point(data = as.data.frame(obs_coords), 
             aes(x = x, y = y), 
             shape = 4, size = 2, color = "red", stroke = 1.5) +
  coord_sf(xlim = xlim_orig, ylim = ylim_orig, expand = FALSE) +
  labs(title = paste("Relative Cover Lost -", max(sim_dates)),
       subtitle = "Red X = observed disease presence") +
  theme_classic(base_family = "Georgia")

# Plot 3: Absolute cover lost (grid extent)
plot3 <- ggplot() +
  geom_sf(data = grid_sf, aes(fill = final_R), color = NA) +
  scale_fill_viridis(name = "Cover lost\n(absolute)", option = "plasma", 
                     na.value = "transparent") +
  geom_spatvector(data = land_vect, fill = "gray90", color = "black") +
  geom_point(data = as.data.frame(obs_coords), 
             aes(x = x, y = y), 
             shape = 4, size = 2, color = "red", stroke = 1.5) +
  coord_sf(xlim = xlim_grid, ylim = ylim_grid, expand = FALSE) +
  labs(title = paste("Absolute Cover Lost (Grid Extent) -", max(sim_dates)),
       subtitle = "Red X = observed disease presence") +
  theme_classic(base_family = "Georgia")

# Plot 4: Relative cover lost (grid extent)
plot4 <- ggplot() +
  geom_sf(data = grid_sf, aes(fill = relative_R), color = NA) +
  scale_fill_viridis(name = "Cover lost\n(proportion)", option = "plasma", 
                     na.value = "transparent", limits = c(0, max_relative_rounded)) +
  geom_spatvector(data = land_vect, fill = "gray90", color = "black") +
  geom_point(data = as.data.frame(obs_coords), 
             aes(x = x, y = y), 
             shape = 4, size = 2, color = "red", stroke = 1.5) +
  coord_sf(xlim = xlim_grid, ylim = ylim_grid, expand = FALSE) +
  labs(title = paste("Relative Cover Lost (Grid Extent) -", max(sim_dates)),
       subtitle = "Red X = observed disease presence") +
  theme_classic(base_family = "Georgia")

print(plot1)
print(plot2)
print(plot3)
print(plot4)

# ggsave(here("output", "model_vs_obs_absolute.png"), 
#        plot1, width = 8, height = 6, dpi = 300, bg = "white")
# ggsave(here("output", "model_vs_obs_relative.png"), 
#        plot2, width = 8, height = 6, dpi = 300, bg = "white")
# ggsave(here("output", "model_vs_obs_absolute_grid.png"), 
#        plot3, width = 8, height = 6, dpi = 300, bg = "white")
# ggsave(here("output", "model_vs_obs_relative_grid.png"), 
#        plot4, width = 8, height = 6, dpi = 300, bg = "white")

################################## Interactive leaflet map ##################################

cat("Creating interactive leaflet map...\n")

library(leaflet)

# Create color palette for absolute cover loss
pal_absolute <- colorNumeric(palette = viridis::plasma(256), 
                             domain = grid_sf$final_R,
                             na.color = "transparent")

# Create observation dataframe for leaflet
obs_df <- data.frame(
  lon = obs_coords[, 1],
  lat = obs_coords[, 2],
  location = obs_2019$location,
  date = obs_2019$date
)

# Create leaflet map
leaflet_map <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addTiles(group = "OpenStreetMap") %>%
  addPolygons(data = grid_sf,
              fillColor = ~pal_absolute(final_R),
              fillOpacity = 0.7,
              color = "black",
              weight = 0.5,
              popup = ~paste0("<b>Grid ID:</b> ", unique_ID, 
                              "<br><b>Absolute Cover Lost:</b> ", round(final_R * 100, 2), "%",
                              "<br><b>Relative Cover Lost:</b> ", round(relative_R * 100, 1), "%")) %>%
  addCircleMarkers(data = obs_df,
                   lng = ~lon, lat = ~lat,
                   radius = 5,
                   color = "red",
                   fillColor = "red",
                   fillOpacity = 0.8,
                   stroke = TRUE,
                   weight = 2,
                   popup = ~paste0("<b>Location:</b> ", location,
                                   "<br><b>Date:</b> ", format(date, "%Y-%m-%d"))) %>%
  addLegend("bottomright", pal = pal_absolute, 
            values = grid_sf$final_R,
            title = "Absolute Cover<br>Lost",
            opacity = 1,
            labFormat = labelFormat(transform = function(x) x * 100, suffix = "%")) %>%
  addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
                   options = layersControlOptions(collapsed = FALSE))

print(leaflet_map)

# htmlwidgets::saveWidget(leaflet_map, 
#                         file = here("output", "model_vs_obs_interactive.html"))

################################## Diagnostic leaflet map for disease threshold ##################################

cat("Creating diagnostic map for disease detection...\n")

# Create a version of grid_sf with disease status at observation times
grid_sf$ever_diseased <- sites_vect$ever_diseased[match(grid_sf$unique_ID, sites_vect$unique_ID)]

# For each observation, get the model disease status at that time
obs_df$model_diseased <- obs_matched$model_diseased
obs_df$site_ever_diseased <- obs_matched$site_ever_diseased
obs_df$nearest_site <- obs_matched$nearest_site
obs_df$distance_m <- obs_matched$distance_m

# Create color palette for disease status
pal_disease <- colorFactor(palette = c("gray80", "orange"), 
                           domain = c(FALSE, TRUE),
                           na.color = "transparent")

# Create marker colors for observations
obs_df$marker_color <- ifelse(obs_df$model_diseased, "green", "red")
obs_df$marker_label <- ifelse(obs_df$model_diseased, 
                              "Model detected disease", 
                              "Model missed disease")

# Summary stats for popup
cat("\n=== DISEASE DETECTION SUMMARY ===\n")
cat("Observations where model detected disease:", sum(obs_matched$model_diseased), "/", nrow(obs_matched), "\n")
cat("Mean distance to nearest site:", round(mean(obs_matched$distance_m)), "m\n")
cat("Sites ever diseased:", sum(sites_vect$ever_diseased), "/", length(sites_vect$ever_diseased), "\n\n")

# Create diagnostic leaflet map
diagnostic_map <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
  addTiles(group = "OpenStreetMap") %>%
  
  # Add grid squares colored by "ever diseased" status
  addPolygons(data = grid_sf,
              fillColor = ~pal_disease(ever_diseased),
              fillOpacity = 0.6,
              color = "black",
              weight = 0.5,
              popup = ~paste0("<b>Grid ID:</b> ", unique_ID, 
                              "<br><b>Ever Diseased:</b> ", ever_diseased,
                              "<br><b>Final Cover Lost:</b> ", round(final_R * 100, 3), "%",
                              "<br><b>Threshold:</b> 0.1%")) %>%
  
  # Add observations with color coding
  addCircleMarkers(data = obs_df,
                   lng = ~lon, lat = ~lat,
                   radius = 6,
                   color = ~marker_color,
                   fillColor = ~marker_color,
                   fillOpacity = 0.9,
                   stroke = TRUE,
                   weight = 2,
                   popup = ~paste0("<b>Location:</b> ", location,
                                   "<br><b>Date:</b> ", format(date, "%Y-%m-%d"),
                                   "<br><b>Nearest Site:</b> ", nearest_site,
                                   "<br><b>Distance:</b> ", round(distance_m), " m",
                                   "<br><b>Status:</b> ", marker_label,
                                   "<br><b>Site Ever Diseased:</b> ", site_ever_diseased)) %>%
  
  addLegend("bottomright", 
            colors = c("gray80", "orange", "green", "red"),
            labels = c("Grid: Never diseased", "Grid: Ever diseased", 
                       "Obs: Model detected", "Obs: Model missed"),
            title = "Disease Status<br>(Threshold: 0.1%)",
            opacity = 1) %>%
  
  addLayersControl(baseGroups = c("Satellite", "OpenStreetMap"),
                   options = layersControlOptions(collapsed = FALSE))

print(diagnostic_map)

# htmlwidgets::saveWidget(diagnostic_map, 
#                         file = here("output", "model_disease_diagnostic.html"))

cat("\n✓ Script complete!\n")