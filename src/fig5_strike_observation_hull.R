 
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(rhdf5)
  library(terra)
  library(sf)
  library(tidyverse)
  library(tidyterra)
  library(patchwork)
  library(extrafont)
  library(ggnewscale)
  
  ################################## PLOTTING OPTIONS ##################################
  
  # Spatial extent filter for observations (xmin, xmax, ymin, ymax)
  obs_extent <- ext(-65.50911, -64.57022, 18.04508, 18.64853)
  
  # Plot bounds
  PLOT_BOUNDS <- list(
    xlim = c(obs_extent[1], obs_extent[2]),
    ylim = c(obs_extent[3], obs_extent[4])
  )
  
  # # Plot bounds (single plot, no subplots)
  # PLOT_BOUNDS <- list(
  #   xlim = c(-65.50911, -64.11830),
  #   ylim = c(18.04508, 18.83152)
  # )
  
  # Font settings
  FONT_FAMILY <- "Georgia"
  
  # Text sizes
  TITLE_SIZE <- 11
  AXIS_TEXT_SIZE <- 7
  LEGEND_TEXT_SIZE <- 8
  
  # Figure dimensions
  FIG_WIDTH <- 6
  FIG_HEIGHT <- 4
  FIG_DPI <- 300
  
  # Model site settings
  SITE_SIZE <- 0.08
  SITE_ALPHA <- 1.0
  # SITE_COLOR_DISEASED <- "#D55E00"  # Orange (colorblind-friendly)
  # SITE_COLOR_DISEASED <- "#E6B800"
  SITE_COLOR_DISEASED <- "#4B0082"
  SITE_COLOR_HEALTHY <- "gray80"
  
  # Observation marker settings
  OBS_SIZE <- 0.7
  OBS_STROKE <- 0.1
  OBS_ALPHA <- 0.6
  OBS_COLOR_PRESENCE <- "red"
  OBS_COLOR_ABSENCE <- "gray30"
  
  # Hull settings
  HULL_COLOR <- "red"
  HULL_LINEWIDTH <- 0.5
  HULL_LINETYPE <- "dashed"
  HULL_ALPHA <- 0.6
  
  # Land settings
  LAND_FILL <- "gray90"
  LAND_COLOR <- "black"
  
  # Background settings
  # OCEAN_BG_COLOR <- "gray80"
  OCEAN_BG_COLOR <- "white"
  
  # Disease threshold (>0.1% cover lost = diseased)
  DISEASE_THRESH <- 0.001
  
  # Legend dot size (override for visibility)
  LEGEND_DOT_SIZE <- 2
  
  ################################## Load Data ##################################
  
  matlab_file <- here("output", "seascape_SIR", "seascape_SIR_workspace.mat")
  
  locations_raw <- h5read(matlab_file, "/outputs/sites/locations")
  locations <- if (ncol(locations_raw) > nrow(locations_raw)) t(locations_raw) else locations_raw
  
  unique_IDs <- as.vector(h5read(matlab_file, "/outputs/sites/unique_IDs"))
  R_total <- h5read(matlab_file, "/outputs/totals/R_total")
  tspan <- as.vector(h5read(matlab_file, "/outputs/metadata/tspan"))
  
  ref_date <- as.Date("2018-12-01")
  sim_dates <- ref_date + (tspan[1]:tspan[2]) - 1
  
  obs <- read.csv(here("output", "combined_coral_data.csv")) %>%
    mutate(date = as.Date(date),
           quarter = case_when(
             date < as.Date("2019-01-01") ~ "Q1",
             date < as.Date("2019-04-01") ~ "Q1",
             date < as.Date("2019-07-01") ~ "Q2",
             date < as.Date("2019-10-01") ~ "Q3",
             date < as.Date("2020-01-01") ~ "Q4",
             TRUE ~ NA_character_
           ),
           obs_type = case_when(
             presence %in% c("P", "S") ~ "Presence",
             presence == "A" ~ "Absence",
             TRUE ~ NA_character_
           )) %>%
    filter(!is.na(quarter), !is.na(obs_type)) %>%
    filter(lon >= obs_extent[1] & lon <= obs_extent[2] &
             lat >= obs_extent[3] & lat <= obs_extent[4])
  
  land_vect <- readRDS(here("output", "osm_land_vect.rds"))
  
  ################################## Load Bathymetry for Background ##################################
  
  spatial_metadata <- readRDS(here('data/spatial_metadata.rds'))
  bathy_final <- readRDS(here('data/bathy_final.rds'))
  terra::crs(bathy_final) <- spatial_metadata$bathy_final$crs
  
  # Create depth mask for 0-60m
  depth_mask <- bathy_final >= -60 & bathy_final <= 0
  depth_zone_raster <- depth_mask
  depth_zone_raster[depth_mask == 0] <- NA
  depth_zone_raster[depth_mask == 1] <- 1
  
  # Project to lat/lon
  depth_zone_latlon <- project(depth_zone_raster, "EPSG:4326", method = "near")
  
  ################################## Create Quarterly Plots ##################################
  
  quarters <- c("Q1", "Q2", "Q3", "Q4")
  quarter_dates <- as.Date(c("2019-03-31", "2019-06-30", "2019-09-30", "2019-12-31"))
  
  plots <- map2(quarters, quarter_dates, function(qtr, end_date) {
    
    day_idx <- max(which(sim_dates <= end_date), 1)
    disease_status <- apply(R_total[1:day_idx, , drop = FALSE], 2, max) >= DISEASE_THRESH
    
    sites_df <- data.frame(lon = locations[, 1], lat = locations[, 2], diseased = disease_status)
    
    obs_qtr <- obs %>% filter(quarter == qtr)
    obs_presence <- obs_qtr %>% filter(obs_type == "Presence")
    
    hull <- if (nrow(obs_presence) >= 3) {
      st_convex_hull(st_union(st_as_sf(obs_presence, coords = c("lon", "lat"), crs = 4326)))
    } else NULL
    
    # Crop depth raster to plot bounds
    buffer <- 0.1
    region_extent <- ext(PLOT_BOUNDS$xlim[1]-buffer, PLOT_BOUNDS$xlim[2]+buffer,
                         PLOT_BOUNDS$ylim[1]-buffer, PLOT_BOUNDS$ylim[2]+buffer)
    depth_cropped <- crop(depth_zone_latlon, region_extent)
    
    ggplot() +
      # Depth zone background (gray25)
      geom_spatraster(data = depth_cropped) +
      scale_fill_gradient(low = OCEAN_BG_COLOR, high = OCEAN_BG_COLOR, 
                          na.value = "transparent", guide = "none") +
      # Model sites
      geom_point(data = sites_df, aes(lon, lat, color = diseased), 
                 size = SITE_SIZE, alpha = SITE_ALPHA) +
      scale_color_manual(values = c("FALSE" = SITE_COLOR_HEALTHY, "TRUE" = SITE_COLOR_DISEASED),
                         labels = c("Unaffected", "SCTLD-affected"), name = NULL) +
      # Land
      geom_spatvector(data = land_vect, fill = LAND_FILL, color = LAND_COLOR, linewidth = 0.3) +
      # Hull
      {if (!is.null(hull)) geom_sf(data = hull, fill = NA, color = HULL_COLOR, 
                                   linewidth = HULL_LINEWIDTH, linetype = HULL_LINETYPE,
                                   alpha = HULL_ALPHA)} +
      # Reset fill scale for observations
      new_scale_fill() +
      # Observations
      geom_point(data = obs_qtr, aes(lon, lat, fill = obs_type), 
                 shape = 21, size = OBS_SIZE, stroke = OBS_STROKE, alpha = OBS_ALPHA) +
      scale_fill_manual(values = c("Presence" = OBS_COLOR_PRESENCE, "Absence" = OBS_COLOR_ABSENCE),
                        name = NULL) +
      coord_sf(xlim = PLOT_BOUNDS$xlim, ylim = PLOT_BOUNDS$ylim, expand = FALSE) +
      labs(x = NULL, y = NULL) +
      theme_classic(base_family = FONT_FAMILY, base_size = AXIS_TEXT_SIZE) +
      theme(legend.position = "none",
            axis.text = element_text(size = AXIS_TEXT_SIZE)) +
      guides(color = guide_legend(override.aes = list(size = LEGEND_DOT_SIZE)),
             fill = guide_legend(override.aes = list(size = LEGEND_DOT_SIZE)))
  })
  
  ################################## Combine and Save ##################################
  
  combined <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom",
          legend.text = element_text(size = LEGEND_TEXT_SIZE, family = FONT_FAMILY),
          plot.tag = element_text(size = TITLE_SIZE, family = FONT_FAMILY),
          plot.tag.position = c(0.95, 0.95))
  
  print(combined)
  
  ggsave(here("output", "output_figures_tables", "quarterly_disease_presence.png"), 
         combined, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI, bg = "white")