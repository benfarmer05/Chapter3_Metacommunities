# .rs.restartR(clean = TRUE)
rm(list=ls())

library(here)
library(rhdf5)
library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(patchwork)

################################## CONTROL TOGGLE ##################################

# Spatial extent filter for observations (xmin, xmax, ymin, ymax)
obs_extent <- ext(-65.58409, -64.50770, 18.10861, 18.64853)

# Plot extent (xmin, xmax, ymin, ymax) - defaults to obs_extent
plot_xlim <- c(obs_extent[1], obs_extent[2])
plot_ylim <- c(obs_extent[3], obs_extent[4])

# Text settings
text_family <- "Georgia"
text_size_title <- 11
text_size_axis <- 9
text_size_legend <- 9

# Point settings
point_size <- 1.2
point_alpha <- 0.7
obs_marker_size <- 2.5
obs_marker_stroke <- 1.2

# Hull settings
hull_color <- "blue"
hull_linewidth <- 1.5
hull_linetype <- "dashed"

# Color settings
diseased_color <- "red"
healthy_color <- "gray80"
obs_color <- "blue"
land_fill <- "gray90"
land_color <- "black"

# Legend settings
legend_position <- "right"

################################################################################

# Load MATLAB outputs
matlab_file <- here("output", "seascape_SIR", "seascape_SIR_workspace.mat")

# Read locations and transpose if needed
locations_raw <- h5read(matlab_file, "/outputs/sites/locations")
if (ncol(locations_raw) > nrow(locations_raw)) {
  locations <- t(locations_raw)  # Now rows = sites, cols = [lon, lat]
} else {
  locations <- locations_raw
}

unique_IDs <- as.vector(h5read(matlab_file, "/outputs/sites/unique_IDs"))
R_total <- h5read(matlab_file, "/outputs/totals/R_total")
tspan <- as.vector(h5read(matlab_file, "/outputs/metadata/tspan"))

# Create date sequence
ref_date <- as.Date("2018-12-01")
sim_dates <- ref_date + (tspan[1]:tspan[2]) - 1

# Disease threshold (>0.1% cover lost = diseased)
disease_thresh <- 0.001

# Load observations
obs <- read.csv(here("output", "combined_coral_data.csv")) %>%
  mutate(date = as.Date(date),
         quarter = case_when(
           date < as.Date("2019-01-01") ~ "Q1",  # Dec 2018 -> Q1
           date < as.Date("2019-04-01") ~ "Q1",
           date < as.Date("2019-07-01") ~ "Q2",
           date < as.Date("2019-10-01") ~ "Q3",
           date < as.Date("2020-01-01") ~ "Q4",  # Q4 2019 only
           TRUE ~ NA_character_  # Everything else (2020+)
         )) %>%
  filter(presence %in% c("P", "S")) %>%
  filter(!is.na(quarter)) %>%  # Keep only 2019 data (Q1-Q4)
  # Filter to desired spatial extent
  filter(lon >= obs_extent[1] & lon <= obs_extent[2] &
           lat >= obs_extent[3] & lat <= obs_extent[4])

# Load land (terra SpatVector)
land_vect <- readRDS(here("output", "osm_land_vect.rds"))

# Calculate quarterly disease status
quarters <- c("Q1", "Q2", "Q3", "Q4")
quarter_dates <- as.Date(c("2019-03-31", "2019-06-30", "2019-09-30", "2019-12-31"))

plots <- map2(quarters, quarter_dates, function(qtr, end_date) {
  # Find last day of quarter
  day_idx <- which(sim_dates <= end_date)
  if (length(day_idx) == 0) day_idx <- 1
  day_idx <- max(day_idx)
  
  # Get cumulative disease status (ever diseased by end of quarter)
  disease_status <- apply(R_total[1:day_idx, , drop = FALSE], 2, max) >= disease_thresh
  
  # Create site dataframe
  sites_df <- data.frame(
    lon = locations[, 1],
    lat = locations[, 2],
    diseased = disease_status
  )
  
  # Filter observations for this quarter
  obs_qtr <- obs %>% filter(quarter == qtr)
  
  # Create convex hull around observations
  if (nrow(obs_qtr) >= 3) {
    obs_sf <- st_as_sf(obs_qtr, coords = c("lon", "lat"), crs = 4326)
    hull <- st_convex_hull(st_union(obs_sf))
  } else {
    hull <- NULL
  }
  
  # Create plot
  p <- ggplot() +
    geom_point(data = sites_df, aes(lon, lat, color = diseased), 
               size = point_size, alpha = point_alpha) +
    scale_color_manual(
      values = c("FALSE" = healthy_color, "TRUE" = diseased_color),
      labels = c("Inactive/Undetectable", "Observable SCTLD loss"),
      name = NULL
    ) +
    geom_spatvector(data = land_vect, fill = land_fill, color = land_color, linewidth = 0.3) +
    {if (!is.null(hull)) geom_sf(data = hull, fill = NA, color = hull_color, 
                                 linewidth = hull_linewidth, linetype = hull_linetype)} +
    geom_point(data = obs_qtr, aes(lon, lat), 
               shape = 4, size = obs_marker_size, color = obs_color, stroke = obs_marker_stroke) +
    coord_sf(xlim = plot_xlim, ylim = plot_ylim, expand = FALSE) +
    labs(title = paste(qtr, "2019"), 
         subtitle = paste(nrow(obs_qtr), "observations"),
         x = NULL, y = NULL) +
    theme_classic(base_family = text_family, base_size = text_size_axis) +
    theme(
      legend.position = "none",  # Will combine legends later
      plot.title = element_text(hjust = 0.5, face = "bold", size = text_size_title),
      plot.subtitle = element_text(hjust = 0.5, size = text_size_axis),
      axis.text = element_text(size = text_size_axis)
    )
  
  return(p)
})

# Combine plots with shared legend on the right
combined <- wrap_plots(plots, ncol = 2, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(
    legend.position = legend_position,
    legend.text = element_text(size = text_size_legend, family = text_family),
    legend.title = element_text(size = text_size_legend, family = text_family, face = "bold")
  )

print(combined)

# Save
# ggsave(here("output", "quarterly_disease_presence.png"), 
#        combined, width = 12, height = 10, dpi = 300, bg = "white")