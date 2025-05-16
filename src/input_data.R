  
  # .rs.restartR(clean = TRUE)
  rm(list=ls())
  
  library(here)
  library(readxl)
  library(tidyverse)
  library(gganimate)
  library(MASS)
  library(sf)
  library(ggspatial)
  
  ################################## Set-up ##################################
  
  # hunt = read_excel(here("data/Hunt Backup Spreadsheet.xlsx"))
  # interv = read_excel(here("data/SCTLD Intervention Information.xlsx"),
  #                     sheet = 'Intervention')
  rover = read.csv(here("data/SCTLDRoverSurveyFULLspecies_0.csv"))
  emerge = read_excel(here("data/SCTLDSites_DateEmerged_Oct2020.xlsx"))
  
  load(here("data/USVI_2019_benthic_cover.rda"))
  
  #clean up dataframes
  rover = rover %>%
    rename(island = Island, team = Team, date = Date, presence = SCTLD.Present.Absent.Suspect...P.A.S., type = Survey.Type..HFCD..Roving.Survey..Researcher.Report.,
           location = Site.Name.or.Landmark, surveyors = Surveyor.s., comorb = Other.sig..impairments.observed,
           severity = Level.of.severity.of.SCTLD, lat = Latitude, lon = Longitude) %>%
    dplyr::select(-Grid..optional., -Zone..optional., -Agg..Spp., -DLAB, -DSTO, -DCYL, -EFAS, -PCLI, -PSTR, -CNAT, -MMEA, -MCAV, -OFAV,
           -OFRA, -OANN, -SSID, -SBOU, -SINT, -Data.Entry..Initials., -Monitoring...Y.N., -Data.sheet.scanned...Initials.,
           -Photos.uploaded...initals., -ObjectId, -x, -y, -Intervention.needed...Y.N., -Type.Recommendation,
           -Date.of.Intervention, -Type.Intervention.done, -surveyors) %>%
    mutate(date = as.POSIXct(date, format = '%m/%d/%Y %I:%M:%S %p', tz = 'UTC')) %>%
    # filter(presence == 'P' | presence == 'S') %>%
    filter(island == 'STT' | island == 'STJ') %>%
    mutate(island = as.factor(island)) %>%
    filter(year(date) == 2019) %>%
    arrange(date)
  
  #remove species-specific qualitative prevalence info if it isn't needed
  rover = rover %>%
    dplyr::select(-MMEA.Prevalence, -DCYL.Prevalence, -DSTO.Prevalence, -EFAS.Prevalence, -CNAT.Prevalence, -DLAB.Prevalence,
           -PSTR.Prevalence, -PCLI.Prevalence)
  
  emerge = emerge %>%
    rename(location = Location, lat = Latitude, lon = Longitude, date = `Date First emerged`) %>%
    filter(grepl("^\\d{4}_[A-Za-z]{3}$", `date`)) %>% #filters out Coral Bay, which hadn't been recorded with any infections at time of this worksheet's construction
    mutate( #adjust date to final day of month first observed as emerged, for conservatism
      date = parse_date_time(`date`, "Y_b"),
      date = as.POSIXct(rollback(date + months(1)))
    ) %>%
    filter(year(date) == 2018 | year(date) == 2019) %>%
    arrange(date)
  
  #update dates associated with NCRMP data in 'rover' to more accurately represent the day in which SCTLD presence was noted during
  #   2019 NCRMP mission
  # Assume 'rover' and 'USVI_2019_benthic_cover' are already in your environment
  # Step 1: Prepare a lookup table of PSU codes to POSIXct dates
  cover_dates <- USVI_2019_benthic_cover %>%
    transmute(
      psu = as.character(`PRIMARY_SAMPLE_UNIT`),        # make sure code is character
      obs_date = make_datetime(year = YEAR,               # build POSIXct from YEAR, MONTH, DAY
                               month = MONTH,
                               day = DAY,
                               tz = "UTC")
    ) %>%
    distinct(psu, obs_date) %>%              # remove exact duplicates
    group_by(psu) %>%                        # if multiple dates exist per PSU
    summarize(obs_date = min(obs_date),      # choose earliest (or change to max() as needed)
              .groups = "drop")            # one row per psu
  
  # Step 2: Extract 4-digit NCRMP code from rover$location, join, and update date
  rover <- rover %>%
    mutate(
      # extract the 4-digit code when location contains "NCRMP ####"
      psu = str_extract(location, "(?<=NCRMP )\\d{4}"),
      # ensure original date is POSIXct for consistent type
      date = as.POSIXct(date)
    ) %>%
    left_join(cover_dates, by = "psu") %>%         # bring in obs_date by matching psu
    mutate(
      # if we have a matched obs_date, replace date; otherwise keep original
      date = if_else(!is.na(obs_date), obs_date, date)
    ) %>%
    dplyr::select(-psu, -obs_date) %>%                        # remove helper columns
    arrange(date)
  
  ################################## Map Emergence spread ##################################

  # Convert to sf
  emerge_sf <- st_as_sf(emerge, coords = c("lon", "lat"), crs = 4326)
  
  # Extract month and year for coloring
  emerge_sf$month <- format(emerge_sf$date, "%b %Y")  # Month-Year as factor
  
  # Convert 'month' to a factor with levels ordered chronologically
  emerge_sf$month <- factor(emerge_sf$month, levels = unique(format(emerge_sf$date, "%b %Y")))
  
  # Bounding box for basemap
  bbox <- st_bbox(emerge_sf)
  
  # Adjusting the bounding box with a larger buffer
  bbox <- st_bbox(emerge_sf)
  buffer_dist <- 0.1  # 0.1 degree buffer around the points
  
  # Increase the bounding box size to avoid small map tiles issue
  bbox_extended <- c(
    xmin = bbox["xmin"] - buffer_dist,
    xmax = bbox["xmax"] + buffer_dist,
    ymin = bbox["ymin"] - buffer_dist,
    ymax = bbox["ymax"] + buffer_dist
  )
  
  # Plot
  ggplot() +
    annotation_map_tile(type = "osm", zoomin = 0) +  # basemap tile
    geom_sf(data = emerge_sf, aes(color = month), size = 3) +
    scale_color_viridis_d() +  # Discrete color scale for months
    theme_minimal() +
    labs(title = "Emerge Locations by Month over OpenStreetMap",
         color = "Month",
         caption = "Basemap: OpenStreetMap") +
    coord_sf(xlim = c(bbox_extended["xmin"], bbox_extended["xmax"]),
             ylim = c(bbox_extended["ymin"], bbox_extended["ymax"]))
  
  ################################## Map Roving spread ##################################
  
  # coral_2019 = rover %>%
  #   filter(year(date) == 2019)
  coral_2019 = rover
  
  # Create a count of occurrences by location for the heatmap intensity
  location_counts <- coral_2019 %>%
    group_by(lat, lon, island) %>%
    summarize(count = n(), first_date = min(date), .groups = 'drop') %>%
    arrange(first_date)
  
  # Basic static map showing all 2019 disease occurrences with size based on counts
  static_map <- ggplot(location_counts, aes(x = lon, y = lat)) +
    borders("world", regions = "Puerto Rico", fill = "gray90", colour = "gray70") +
    geom_point(aes(size = count, color = island, alpha = 0.7)) +
    scale_color_viridis_d() +
    scale_size_continuous(range = c(3, 10)) +
    coord_fixed(xlim = c(-65.2, -64.5), ylim = c(17.6, 18.5)) +
    labs(title = "Coral Disease Occurrences in USVI (2019)",
         subtitle = "Size indicates number of reported cases",
         x = "Longitude", y = "Latitude") +
    theme_minimal()
  
  #assign labels by month
  coral_2019$month <- month(coral_2019$date)
  coral_2019$month_label <- format(coral_2019$date, "%b %Y")
  
  # Create animated map showing progression over time
  animated_map <- ggplot(coral_2019, aes(x = lon, y = lat)) +
    borders("world", regions = "Puerto Rico", fill = "gray90", colour = "gray70") +
    geom_point(aes(color = island, size = 5), alpha = 0.7) +
    scale_color_viridis_d() +
    coord_fixed(xlim = c(-65.2, -64.5), ylim = c(17.6, 18.5)) +
    labs(title = "Coral Disease Spread in USVI",
         subtitle = "Date: {frame_time}",
         x = "Longitude", y = "Latitude") +
    theme_minimal() +
    transition_time(date) +
    shadow_mark(past = TRUE, alpha = 0.3)
  
  
  # coral_2019$month_label <- factor(coral_2019$month_label, 
  #                                  levels = coral_2019$month_label[order(coral_2019$month)])
  
  # coral_2019$month_label <- factor(coral_2019$month_label)
  
  coral_2019$month <- as.numeric(coral_2019$month)
  coral_2019 <- coral_2019[order(coral_2019$month), ]
  coral_2019$month_label <- factor(coral_2019$month_label, 
                                   levels = unique(coral_2019$month_label))
  
  # Get only the months that are actually present
  present_months <- sort(unique(coral_2019$month))
  present_labels <- levels(coral_2019$month_label)
  
  
  # Create a static map showing all disease occurrences
  static_map <- ggplot(coral_2019, aes(x = lon, y = lat)) +
    # Base map - you might need to adjust borders based on your exact location
    borders("world", regions = c("Puerto Rico", "Virgin Islands, U.S."), 
            fill = "gray90", colour = "gray70") +
    
    # Plot points with color by month and shape by severity
    geom_point(aes(color = month, shape = severity), size = 3, alpha = 0.7) + #color = island, 
    
    scale_color_gradientn(
      colours = c("blue", "cyan", "yellow", "orange", "red"),
      breaks = present_months,
      labels = present_labels,
      name = "Month"
    ) +
    
    # Set map boundaries based on your data
    coord_fixed(xlim = c(min(coral_2019$lon) - 0.05, max(coral_2019$lon) + 0.05), 
                ylim = c(min(coral_2019$lat) - 0.05, max(coral_2019$lat) + 0.05)) +
    
    # Labels
    labs(
      title = "Coral Disease Occurrences in USVI (2019)",
      subtitle = "All reported cases",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal()
  
  # # Save the static map
  # ggsave(here("output/coral_disease_all_occurrences_2019.png"), static_map, width = 10, height = 8)
  
  
  
  

  # Add month for temporal grouping
  coral_2019 <- coral_2019 %>%
    arrange(date) %>%
    mutate(
      month = floor_date(date, "month"),
      month_label = format(month, "%b %Y")
    )
  
  # Create monthly cumulative datasets for animation and facets
  # Get unique months in the dataset
  unique_months <- unique(coral_2019$month)
  
  # Function to create a safe kernel density estimation that won't fail
  # even with small numbers of points
  safe_kde <- function(points_data, bw = 0.01) {
    if(nrow(points_data) < 5) {
      # Not enough points for reliable KDE, return a simple grid
      x_range <- seq(min(coral_2019$lon) - 0.05, max(coral_2019$lon) + 0.05, length.out = 100)
      y_range <- seq(min(coral_2019$lat) - 0.05, max(coral_2019$lat) + 0.05, length.out = 100)
      grid <- expand.grid(x = x_range, y = y_range)
      grid$z <- 0
      return(grid)
    } else {
      # Enough points for KDE
      kde <- kde2d(
        points_data$lon, 
        points_data$lat,
        n = 100,
        lims = c(min(coral_2019$lon) - 0.05, max(coral_2019$lon) + 0.05, 
                 min(coral_2019$lat) - 0.05, max(coral_2019$lat) + 0.05),
        h = bw
      )
      
      grid <- expand.grid(x = kde$x, y = kde$y)
      grid$z <- as.vector(kde$z)
      return(grid)
    }
  }
  
  # Create cumulative datasets for each month
  monthly_kde_list <- list()
  monthly_points_list <- list()
  
  for(i in 1:length(unique_months)) {
    current_month <- unique_months[i]
    
    # Get all points up to and including this month
    cumulative_points <- coral_2019 %>%
      filter(date <= ceiling_date(current_month, "month") - days(1))
    
    # Store points data for this month
    monthly_points_list[[i]] <- cumulative_points %>%
      mutate(month_idx = i, month_date = current_month)
    
    # Generate KDE for these points
    kde_grid <- safe_kde(cumulative_points)
    kde_grid$month_idx <- i
    kde_grid$month_date <- current_month
    kde_grid$month_label <- format(current_month, "%b %Y")
    
    monthly_kde_list[[i]] <- kde_grid
  }
  
  # Combine all KDE grids
  all_kde <- do.call(rbind, monthly_kde_list)
  
  # Combine all points
  all_points <- do.call(rbind, monthly_points_list)
  
  
  
  #test
  ordered_labels <- format(sort(unique_months), "%b %Y")
  all_kde$month_label <- factor(all_kde$month_label, levels = ordered_labels)
  all_points$month_label <- factor(format(all_points$month_date, "%b %Y"), levels = ordered_labels)
  
  
  # Create a faceted map showing cumulative spread by month
  faceted_map <- ggplot() +
    # Add KDE layer
    geom_tile(data = all_kde, aes(x = x, y = y, fill = z), alpha = 0.7) +
    
    # Add points layer
    geom_point(data = all_points, aes(x = lon, y = lat, color = island), size = 1.5, alpha = 0.7) +
    
    # Colors
    scale_fill_viridis_c(option = "plasma", name = "Density") +
    scale_color_brewer(palette = "Set1") +
    
    # Facet by month
    facet_wrap(~month_label, ncol = 3) +
    
    # Map bounds - adjust as needed for your data
    coord_fixed() +
    
    # Labels
    labs(
      title = "Cumulative Coral Disease Spread in USVI (2019)",
      subtitle = "Each panel shows cases accumulated through that month",
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "steelblue", color = "steelblue"), #, alpha = 0.8
      strip.text = element_text(color = "white", face = "bold")
    )
  
  # # Save the faceted map
  # ggsave(here("output/coral_disease_monthly_facets_2019.png"), faceted_map, width = 12, height = 10)
  

  ################################## Combined spread ##################################
  
  # Add a source column to each dataframe
  emerge <- emerge %>%
    mutate(source = "emerge")
  
  rover <- rover %>%
    mutate(source = "rover")
  
  # Combine the two dataframes
  combined <- bind_rows(emerge %>% dplyr::select(location, lat, lon, date, source),
                        rover %>% dplyr::select(location, lat, lon, date, source)) %>%
    arrange(date)
  