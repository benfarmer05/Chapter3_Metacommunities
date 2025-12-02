  
  # rm(list = ls())
  
  library(here)
  library(rhdf5)
  library(tidyverse)
  library(patchwork)
  library(cowplot)
  library(extrafont)
  
  # font_import()  # Only needed once - imports all system fonts
  # loadfonts()    # Load fonts for current session
  
  ################################## Settings ##################################
  
  use_existing_cache <- FALSE  # Set to TRUE to load cached data, FALSE to regenerate
  use_group_denominator <- FALSE  # TRUE = divide by group total, FALSE = divide by site total
  
  # Plot styling parameters
  linewidths <- 0.2   # Line width in mm (0.75 mm ≈ 1 pt)
  titlesize <- 10     # Title text size in points
  textsize <- 9       # Axis text size in points
  linealpha <- 0.4    # Transparency for individual site lines
  legend_linewidth <- 1.5  # Line width for legend lines (in mm)
  
  # Figure dimensions (in inches)
  fig_width_scen1 <- 7.087
  fig_height_scen1 <- 5
  fig_width_scen2_3 <- 7.087
  fig_height_scen2_3 <- 5
  
  # Susceptibility group colors (matching reference)
  group_colors <- c(
    "LS" = "#1E90FF",   # Blue (Low Susceptibility)
    "MS" = "#FFD700",   # Yellow (Moderate Susceptibility)
    "HS" = "#FF1493"    # Deep Pink (High Susceptibility)
  )
  
  seascape_dir <- here("output", "seascape_SIR")
  all_folders <- list.dirs(seascape_dir, full.names = FALSE, recursive = FALSE)
  scenarios <- all_folders[grepl("SCENARIO", all_folders, ignore.case = TRUE)]
  scenarios <- sort(scenarios)
  
  cat("Detected", length(scenarios), "scenario folders:\n")
  for (s in scenarios) cat("  -", s, "\n")
  cat("\n")
  
  if (!dir.exists(here("temp"))) dir.create(here("temp"))
  
  ################################## Data Preparation / Caching ##################################
  
  scenario_data <- list()
  
  for (scenario in scenarios) {
    cat("Processing scenario:", scenario, "\n")
    
    scen_num <- str_extract(scenario, "(?<=SCENARIO)\\d+")
    cache_file <- here("temp", paste0("sir_data_", gsub("[^A-Za-z0-9]", "_", scenario), ".rds"))
    
    if (use_existing_cache && file.exists(cache_file)) {
      cat("  Loading from cache...\n")
      scenario_data[[scenario]] <- readRDS(cache_file)
      
    } else {
      matlab_file <- here("output", "seascape_SIR", scenario, "seascape_SIR_workspace.mat")
      
      if (!file.exists(matlab_file)) {
        cat("  WARNING: File not found, skipping\n")
        next
      }
      
      # Load all compartments
      S_total <- h5read(matlab_file, "/outputs/totals/S_total")
      I_total <- h5read(matlab_file, "/outputs/totals/I_total")
      R_total <- h5read(matlab_file, "/outputs/totals/R_total")
      S_LS <- h5read(matlab_file, "/outputs/SIR/S_LS")
      S_MS <- h5read(matlab_file, "/outputs/SIR/S_MS")
      S_HS <- h5read(matlab_file, "/outputs/SIR/S_HS")
      I_LS <- h5read(matlab_file, "/outputs/SIR/I_LS")
      I_MS <- h5read(matlab_file, "/outputs/SIR/I_MS")
      I_HS <- h5read(matlab_file, "/outputs/SIR/I_HS")
      R_LS <- h5read(matlab_file, "/outputs/SIR/R_LS")
      R_MS <- h5read(matlab_file, "/outputs/SIR/R_MS")
      R_HS <- h5read(matlab_file, "/outputs/SIR/R_HS")
      
      ref_date <- as.Date("2018-12-01")
      tspan <- as.vector(h5read(matlab_file, "/outputs/metadata/tspan"))
      sim_dates <- ref_date + (tspan[1]:tspan[2]) - 1
      
      # Filter sites: keep only those with max R > 0.00001
      max_R_per_site <- apply(R_total, 2, max, na.rm = TRUE)
      keep_sites <- max_R_per_site > 0.00001
      
      cat("  Keeping", sum(keep_sites), "of", length(keep_sites), 
          "sites (removed", sum(!keep_sites), "with max R ≤ 0.00001)\n")
      
      # Apply filter to all datasets
      S_total <- S_total[, keep_sites]
      I_total <- I_total[, keep_sites]
      R_total <- R_total[, keep_sites]
      S_LS <- S_LS[, keep_sites]
      S_MS <- S_MS[, keep_sites]
      S_HS <- S_HS[, keep_sites]
      I_LS <- I_LS[, keep_sites]
      I_MS <- I_MS[, keep_sites]
      I_HS <- I_HS[, keep_sites]
      R_LS <- R_LS[, keep_sites]
      R_MS <- R_MS[, keep_sites]
      R_HS <- R_HS[, keep_sites]
      
      # Convert totals to long format
      I_long <- I_total %>%
        as.data.frame() %>%
        mutate(date = sim_dates) %>%
        pivot_longer(cols = -date, names_to = "site", values_to = "infected") %>%
        mutate(site = as.numeric(gsub("V", "", site)))
      
      R_long <- R_total %>%
        as.data.frame() %>%
        mutate(date = sim_dates) %>%
        pivot_longer(cols = -date, names_to = "site", values_to = "removed") %>%
        mutate(site = as.numeric(gsub("V", "", site)))
      
      # Convert susceptibility groups to long format
      I_by_group <- bind_rows(
        I_LS %>% as.data.frame() %>% mutate(date = sim_dates, group = "LS") %>%
          pivot_longer(cols = -c(date, group), names_to = "site", values_to = "infected"),
        I_MS %>% as.data.frame() %>% mutate(date = sim_dates, group = "MS") %>%
          pivot_longer(cols = -c(date, group), names_to = "site", values_to = "infected"),
        I_HS %>% as.data.frame() %>% mutate(date = sim_dates, group = "HS") %>%
          pivot_longer(cols = -c(date, group), names_to = "site", values_to = "infected")
      ) %>%
        mutate(site = as.numeric(gsub("V", "", site)))
      
      R_by_group <- bind_rows(
        R_LS %>% as.data.frame() %>% mutate(date = sim_dates, group = "LS") %>%
          pivot_longer(cols = -c(date, group), names_to = "site", values_to = "removed"),
        R_MS %>% as.data.frame() %>% mutate(date = sim_dates, group = "MS") %>%
          pivot_longer(cols = -c(date, group), names_to = "site", values_to = "removed"),
        R_HS %>% as.data.frame() %>% mutate(date = sim_dates, group = "HS") %>%
          pivot_longer(cols = -c(date, group), names_to = "site", values_to = "removed")
      ) %>%
        mutate(site = as.numeric(gsub("V", "", site)))
      
      # Store in list with matrices for relative calculations
      scenario_data[[scenario]] <- list(
        I_long = I_long, 
        R_long = R_long, 
        I_by_group = I_by_group, 
        R_by_group = R_by_group,
        scen_num = scen_num,
        S_total = S_total,
        I_total = I_total,
        R_total = R_total,
        S_LS = S_LS,
        S_MS = S_MS,
        S_HS = S_HS,
        I_LS = I_LS,
        I_MS = I_MS,
        I_HS = I_HS,
        R_LS = R_LS,
        R_MS = R_MS,
        R_HS = R_HS
      )
      
      # Save to cache
      cat("  Saving to cache\n")
      saveRDS(scenario_data[[scenario]], cache_file)
    }
  }
  
  cat("\n✓ Data preparation complete!\n\n")
  
  ################################## Plotting Functions ##################################
  
  plot_totals <- function(data, scen_label) {
    I_long <- data$I_long %>% mutate(infected = infected * 100)
    R_long <- data$R_long %>% mutate(removed = removed * 100)
    
    p_infected <- ggplot(I_long, aes(x = date, y = infected, group = site)) +
      geom_line(alpha = linealpha, color = "black", linewidth = linewidths) +
      scale_y_continuous(labels = scales::comma) +
      scale_x_date(date_labels = "%b-%Y") +
      labs(x = NULL, y = "Coral cover (%)") +
      theme_classic(base_family = "Georgia") +
      theme(
        axis.title.y = element_text(size = titlesize, color = 'black', margin = margin(r = 10)),
        axis.title.x = element_text(size = titlesize, color = 'black'),
        axis.text = element_text(size = textsize, color = 'black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        plot.margin = margin(t = 5, r = 10, b = 0, l = 5),
        plot.tag = element_text(size = titlesize, family = "Georgia"),
        plot.tag.position = c(0.98, 0.98)
      )
    
    p_removed <- ggplot(R_long, aes(x = date, y = removed, group = site)) +
      geom_line(alpha = linealpha, color = "black", linewidth = linewidths) +
      scale_y_continuous(labels = scales::comma) +
      scale_x_date(date_labels = "%b-%Y") +
      labs(x = NULL, y = "Coral cover (%)") +
      theme_classic(base_family = "Georgia") +
      theme(
        axis.title.y = element_text(size = titlesize, color = 'black', margin = margin(r = 10)),
        axis.title.x = element_text(size = titlesize, color = 'black'),
        axis.text = element_text(size = textsize, color = 'black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        plot.margin = margin(t = 5, r = 10, b = 0, l = 5),
        plot.tag = element_text(size = titlesize, family = "Georgia"),
        plot.tag.position = c(0.98, 0.98)
      )
    
    combined_plot <- (p_infected / p_removed) + 
      plot_layout(heights = c(1, 1), axes = "collect_y") +
      plot_annotation(tag_levels = "A")
    
    return(combined_plot)
  }
  
  plot_totals_relative <- function(data, scen_label) {
    S_total <- data$S_total
    I_total <- data$I_total
    R_total <- data$R_total
    
    initial_totals <- data.frame(
      site = 1:ncol(S_total),
      total_init = S_total[1, ] + I_total[1, ] + R_total[1, ]
    )
    
    I_long_rel <- data$I_long %>%
      left_join(initial_totals, by = "site") %>%
      mutate(infected_rel = (infected / total_init) * 100)
    
    R_long_rel <- data$R_long %>%
      left_join(initial_totals, by = "site") %>%
      mutate(removed_rel = (removed / total_init) * 100)
    
    p_infected <- ggplot(I_long_rel, aes(x = date, y = infected_rel, group = site)) +
      geom_line(alpha = linealpha, color = "black", linewidth = linewidths) +
      scale_y_continuous(labels = scales::comma) +
      scale_x_date(date_labels = "%b-%Y") +
      labs(x = NULL, y = "Relative coral cover (%)") +
      theme_classic(base_family = "Georgia") +
      theme(
        axis.title.y = element_text(size = titlesize, color = 'black', margin = margin(r = 10)),
        axis.title.x = element_text(size = titlesize, color = 'black'),
        axis.text = element_text(size = textsize, color = 'black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        plot.margin = margin(t = 5, r = 10, b = 0, l = 5),
        plot.tag = element_text(size = titlesize, family = "Georgia"),
        plot.tag.position = c(0.98, 0.98)
      )
    
    p_removed <- ggplot(R_long_rel, aes(x = date, y = removed_rel, group = site)) +
      geom_line(alpha = linealpha, color = "black", linewidth = linewidths) +
      scale_y_continuous(labels = scales::comma) +
      scale_x_date(date_labels = "%b-%Y") +
      labs(x = NULL, y = "Relative coral cover (%)") +
      theme_classic(base_family = "Georgia") +
      theme(
        axis.title.y = element_text(size = titlesize, color = 'black', margin = margin(r = 10)),
        axis.title.x = element_text(size = titlesize, color = 'black'),
        axis.text = element_text(size = textsize, color = 'black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        plot.margin = margin(t = 5, r = 10, b = 0, l = 5),
        plot.tag = element_text(size = titlesize, family = "Georgia"),
        plot.tag.position = c(0.98, 0.98)
      )
    
    combined_plot <- (p_infected / p_removed) + 
      plot_layout(heights = c(1, 1), axes = "collect_y") +
      plot_annotation(tag_levels = "A")
    
    return(combined_plot)
  }
  
  plot_by_group <- function(data, scen_label) {
    I_by_group <- data$I_by_group %>% 
      mutate(infected = infected * 100,
             group = factor(group, levels = c("LS", "MS", "HS")))
    R_by_group <- data$R_by_group %>% 
      mutate(removed = removed * 100,
             group = factor(group, levels = c("LS", "MS", "HS")))
    
    p_infected_group <- ggplot(I_by_group, aes(x = date, y = infected, 
                                               group = interaction(site, group), 
                                               color = group)) +
      geom_line(alpha = linealpha, linewidth = linewidths) +
      scale_color_manual(
        values = group_colors, 
        labels = c("LS" = "LS", "MS" = "MS", "HS" = "HS"),
        name = "Susc."
      ) +
      scale_y_continuous(labels = scales::comma) +
      scale_x_date(date_labels = "%b-%Y") +
      labs(x = NULL, y = "Coral cover (%)") +
      annotate("text", x = Inf, y = Inf, label = "A", hjust = 1.5, vjust = 1.5,
               size = titlesize / .pt, family = "Georgia") +
      theme_classic(base_family = "Georgia") +
      theme(
        axis.title.y = element_text(size = titlesize, color = 'black', margin = margin(r = 10)),
        axis.title.x = element_text(size = titlesize, color = 'black'),
        axis.text = element_text(size = textsize, color = 'black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = textsize),
        legend.title = element_text(size = titlesize),
        legend.key.height = unit(0.6, "cm"),
        legend.box.spacing = unit(0.3, "cm"),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      ) +
      guides(color = guide_legend(override.aes = list(linewidth = legend_linewidth, alpha = 1)))
    
    p_removed_group <- ggplot(R_by_group, aes(x = date, y = removed, 
                                              group = interaction(site, group), 
                                              color = group)) +
      geom_line(alpha = linealpha, linewidth = linewidths) +
      scale_color_manual(
        values = group_colors,
        labels = c("LS" = "LS", "MS" = "MS", "HS" = "HS"),
        name = "Susc."
      ) +
      scale_y_continuous(labels = scales::comma) +
      scale_x_date(date_labels = "%b-%Y") +
      labs(x = NULL, y = "Coral cover (%)") +
      annotate("text", x = Inf, y = Inf, label = "B", hjust = 1.5, vjust = 1.5,
               size = titlesize / .pt, family = "Georgia") +
      theme_classic(base_family = "Georgia") +
      theme(
        axis.title.y = element_text(size = titlesize, color = 'black', margin = margin(r = 10)),
        axis.title.x = element_text(size = titlesize, color = 'black'),
        axis.text = element_text(size = textsize, color = 'black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = textsize),
        legend.title = element_text(size = titlesize),
        legend.key.height = unit(0.6, "cm"),
        legend.box.spacing = unit(0.3, "cm"),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      ) +
      guides(color = guide_legend(override.aes = list(linewidth = legend_linewidth, alpha = 1)))
    
    combined_plot <- (p_infected_group / p_removed_group) +
      plot_layout(heights = c(1, 1), axes = "collect_y", guides = "collect")
    
    return(combined_plot)
  }
  
  plot_by_group_relative <- function(data, scen_label) {
    S_total <- data$S_total
    I_total <- data$I_total
    R_total <- data$R_total
    S_LS <- data$S_LS
    S_MS <- data$S_MS
    S_HS <- data$S_HS
    I_LS <- data$I_LS
    I_MS <- data$I_MS
    I_HS <- data$I_HS
    R_LS <- data$R_LS
    R_MS <- data$R_MS
    R_HS <- data$R_HS
    
    # Calculate initial totals based on toggle
    if (use_group_denominator) {
      # Divide by group-specific initial total
      initial_totals <- bind_rows(
        data.frame(site = 1:ncol(S_LS), group = "LS", total_init = S_LS[1, ] + I_LS[1, ] + R_LS[1, ]),
        data.frame(site = 1:ncol(S_MS), group = "MS", total_init = S_MS[1, ] + I_MS[1, ] + R_MS[1, ]),
        data.frame(site = 1:ncol(S_HS), group = "HS", total_init = S_HS[1, ] + I_HS[1, ] + R_HS[1, ])
      )
    } else {
      # Divide by site total initial cover (same for all groups at a site)
      site_totals <- data.frame(
        site = 1:ncol(S_total),
        total_init = S_total[1, ] + I_total[1, ] + R_total[1, ]
      )
      initial_totals <- bind_rows(
        site_totals %>% mutate(group = "LS"),
        site_totals %>% mutate(group = "MS"),
        site_totals %>% mutate(group = "HS")
      )
    }
    
    I_by_group_rel <- data$I_by_group %>%
      left_join(initial_totals, by = c("site", "group")) %>%
      mutate(infected_rel = (infected / total_init) * 100,
             group = factor(group, levels = c("LS", "MS", "HS")))
    
    R_by_group_rel <- data$R_by_group %>%
      left_join(initial_totals, by = c("site", "group")) %>%
      mutate(removed_rel = (removed / total_init) * 100,
             group = factor(group, levels = c("LS", "MS", "HS")))
    
    p_infected_group <- ggplot(I_by_group_rel, aes(x = date, y = infected_rel, 
                                                   group = interaction(site, group), 
                                                   color = group)) +
      geom_line(alpha = linealpha, linewidth = linewidths) +
      scale_color_manual(
        values = group_colors, 
        labels = c("LS" = "LS", "MS" = "MS", "HS" = "HS"),
        name = "Susc."
      ) +
      scale_y_continuous(labels = scales::comma) +
      scale_x_date(date_labels = "%b-%Y") +
      labs(x = NULL, y = "Relative coral cover (%)") +
      annotate("text", x = Inf, y = Inf, label = "A", hjust = 1.5, vjust = 1.5,
               size = titlesize / .pt, family = "Georgia") +
      theme_classic(base_family = "Georgia") +
      theme(
        axis.title.y = element_text(size = titlesize, color = 'black', margin = margin(r = 10)),
        axis.title.x = element_text(size = titlesize, color = 'black'),
        axis.text = element_text(size = textsize, color = 'black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = textsize),
        legend.title = element_text(size = titlesize),
        legend.key.height = unit(0.6, "cm"),
        legend.box.spacing = unit(0.3, "cm"),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      ) +
      guides(color = guide_legend(override.aes = list(linewidth = legend_linewidth, alpha = 1)))
    
    p_removed_group <- ggplot(R_by_group_rel, aes(x = date, y = removed_rel, 
                                                  group = interaction(site, group), 
                                                  color = group)) +
      geom_line(alpha = linealpha, linewidth = linewidths) +
      scale_color_manual(
        values = group_colors,
        labels = c("LS" = "LS", "MS" = "MS", "HS" = "HS"),
        name = "Susc."
      ) +
      scale_y_continuous(labels = scales::comma) +
      scale_x_date(date_labels = "%b-%Y") +
      labs(x = NULL, y = "Relative coral cover (%)") +
      annotate("text", x = Inf, y = Inf, label = "B", hjust = 1.5, vjust = 1.5,
               size = titlesize / .pt, family = "Georgia") +
      theme_classic(base_family = "Georgia") +
      theme(
        axis.title.y = element_text(size = titlesize, color = 'black', margin = margin(r = 10)),
        axis.title.x = element_text(size = titlesize, color = 'black'),
        axis.text = element_text(size = textsize, color = 'black'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        legend.position = "right",
        legend.text = element_text(size = textsize),
        legend.title = element_text(size = titlesize),
        legend.key.height = unit(0.6, "cm"),
        legend.box.spacing = unit(0.3, "cm"),
        plot.margin = margin(t = 5, r = 5, b = 0, l = 5)
      ) +
      guides(color = guide_legend(override.aes = list(linewidth = legend_linewidth, alpha = 1)))
    
    combined_plot <- (p_infected_group / p_removed_group) +
      plot_layout(heights = c(1, 1), axes = "collect_y", guides = "collect")
    
    return(combined_plot)
  }
  
  ################################## Create Plots ##################################
  
  output_dir <- here("output", "output_figures_tables")
  
  # Plot totals for each scenario
  fig_total_scen1 <- plot_totals(scenario_data[[1]], paste0("Scen. ", scenario_data[[1]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_totals_scenario", scenario_data[[1]]$scen_num, ".png")),
    plot = fig_total_scen1,
    width = fig_width_scen1,
    height = fig_height_scen1,
    dpi = 300,
    bg = "white"
  )
  
  fig_total_scen2 <- plot_totals(scenario_data[[2]], paste0("Scen. ", scenario_data[[2]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_totals_scenario", scenario_data[[2]]$scen_num, ".png")),
    plot = fig_total_scen2,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  fig_total_scen3 <- plot_totals(scenario_data[[3]], paste0("Scen. ", scenario_data[[3]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_totals_scenario", scenario_data[[3]]$scen_num, ".png")),
    plot = fig_total_scen3,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  # Plot by susceptibility group for each scenario
  fig_group_scen1 <- plot_by_group(scenario_data[[1]], paste0("Scen. ", scenario_data[[1]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_by_group_scenario", scenario_data[[1]]$scen_num, ".png")),
    plot = fig_group_scen1,
    width = fig_width_scen1,
    height = fig_height_scen1,
    dpi = 300,
    bg = "white"
  )
  
  fig_group_scen2 <- plot_by_group(scenario_data[[2]], paste0("Scen. ", scenario_data[[2]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_by_group_scenario", scenario_data[[2]]$scen_num, ".png")),
    plot = fig_group_scen2,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  fig_group_scen3 <- plot_by_group(scenario_data[[3]], paste0("Scen. ", scenario_data[[3]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_by_group_scenario", scenario_data[[3]]$scen_num, ".png")),
    plot = fig_group_scen3,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  # Plot relative totals for each scenario
  fig_total_rel_scen1 <- plot_totals_relative(scenario_data[[1]], paste0("Scen. ", scenario_data[[1]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_totals_relative_scenario", scenario_data[[1]]$scen_num, ".png")),
    plot = fig_total_rel_scen1,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  fig_total_rel_scen2 <- plot_totals_relative(scenario_data[[2]], paste0("Scen. ", scenario_data[[2]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_totals_relative_scenario", scenario_data[[2]]$scen_num, ".png")),
    plot = fig_total_rel_scen2,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  fig_total_rel_scen3 <- plot_totals_relative(scenario_data[[3]], paste0("Scen. ", scenario_data[[3]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_totals_relative_scenario", scenario_data[[3]]$scen_num, ".png")),
    plot = fig_total_rel_scen3,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  # Plot relative by susceptibility group for each scenario
  fig_group_rel_scen1 <- plot_by_group_relative(scenario_data[[1]], paste0("Scen. ", scenario_data[[1]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_by_group_relative_scenario", scenario_data[[1]]$scen_num, ".png")),
    plot = fig_group_rel_scen1,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  fig_group_rel_scen2 <- plot_by_group_relative(scenario_data[[2]], paste0("Scen. ", scenario_data[[2]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_by_group_relative_scenario", scenario_data[[2]]$scen_num, ".png")),
    plot = fig_group_rel_scen2,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  fig_group_rel_scen3 <- plot_by_group_relative(scenario_data[[3]], paste0("Scen. ", scenario_data[[3]]$scen_num))
  ggsave(
    filename = file.path(output_dir, paste0("SIR_by_group_relative_scenario", scenario_data[[3]]$scen_num, ".png")),
    plot = fig_group_rel_scen3,
    width = fig_width_scen2_3,
    height = fig_height_scen2_3,
    dpi = 300,
    bg = "white"
  )
  
  cat("\n✓ All figures complete and saved to:", output_dir, "\n")
  
  
  ################################## quick summary for manuscript, of #/infected sites by scenario ##################################
  
  # --- Total infected sites at the end of each scenario ---
  
  get_final_infected_sites <- function(scenario_list) {
    
    results <- map_df(seq_along(scenario_list), function(i) {
      scen_name <- names(scenario_list)[i]
      I_mat <- scenario_list[[i]]$I_total   # matrix: time × sites
      
      final_infected <- I_mat[nrow(I_mat), ]        # last time step
      infected_sites <- sum(final_infected > 0)     # count sites with infection
      
      tibble(
        scenario = scen_name,
        infected_sites_end = infected_sites
      )
    })
    
    return(results)
  }
  
  # Run it
  final_infected_summary <- get_final_infected_sites(scenario_data)
  print(final_infected_summary)
  