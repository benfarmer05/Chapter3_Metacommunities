  
  # rm(list = ls())
  
  library(here)
  library(ggplot2)
  library(patchwork)
  library(extrafont)
  
  # loadfonts()  # Load fonts for current session
  
  ################################## Settings ##################################
  
  # Plot styling parameters
  linewidth <- 0.75      # Line width in mm
  titlesize <- 10        # Title text size in points
  textsize <- 9          # Axis text size in points
  legend_textsize <- 8   # Legend text size in points
  
  # Figure dimensions (in inches)
  fig_width <- 7
  fig_height <- 3
  
  # Scenario colors (colorblind-friendly palette)
  scenario_colors <- c(
    "Scenario 1" = "#E69F00",  # Orange
    "Scenario 2" = "#0072B2",  # Blue
    "Scenario 3" = "#009E73"   # Bluish green
  )
  
  ################################## Panel A: Internal transmission (I0 threshold) ##################################
  
  # Three scenarios with I0 and corresponding tau = I0/10
  scenarios_internal <- data.frame(
    name = c("Scenario 1", "Scenario 2", "Scenario 3"),
    I0 = c(0.0000010455, 0.00000079, 0.000000652),
    tau = c(0.0000010455/10, 0.00000079/10, 0.000000652/10)
  )
  
  # Create sequence of infection fractions (P_i/N_i)
  max_I0 <- max(scenarios_internal$I0)
  frac_infected <- seq(0, max_I0 * 2, length.out = 1000)
  
  # Calculate internal transmission multiplier for each scenario
  results_internal <- data.frame()
  
  for (i in 1:nrow(scenarios_internal)) {
    I0_val <- scenarios_internal$I0[i]
    tau_val <- scenarios_internal$tau[i]
    
    # Calculate the sigmoid multiplier component: 0.5 * (1 + tanh((frac - I0) / tau))
    multiplier <- 0.5 * (1 + tanh((frac_infected - I0_val) / tau_val))
    
    results_internal <- rbind(results_internal, data.frame(
      frac_infected = frac_infected,
      multiplier = multiplier,
      scenario = scenarios_internal$name[i],
      I0 = I0_val
    ))
  }
  
  results_internal$scenario <- factor(results_internal$scenario, levels = scenarios_internal$name)
  
  # Create Panel A
  panel_A <- ggplot(results_internal, aes(x = frac_infected, y = multiplier, color = scenario)) +
    geom_line(linewidth = linewidth) +
    
    # Add vertical lines at I0 values
    geom_vline(data = scenarios_internal, aes(xintercept = I0, color = name), 
               linetype = "dashed", linewidth = linewidth * 0.75, alpha = 0.7) +
    
    # Add horizontal line at 0.5 (inflection point)
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray40", 
               linewidth = linewidth * 0.75) +
    
    # Scales and labels
    scale_x_continuous(
      labels = function(x) {
        # Format with 1 decimal place in scientific notation
        formatted <- sprintf("%.1e", x)
        return(formatted)
      },
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_color_manual(values = scenario_colors) +
    
    labs(
      x = "Fraction of coral cover infected",
      y = expression(paste(italic(ψ), " step function"))
    ) +
    
    # Theme
    theme_classic(base_family = "Georgia") +
    theme(
      axis.title = element_text(size = titlesize, color = 'black'),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.text = element_text(size = textsize, color = 'black'),
      axis.ticks = element_line(color = "black"),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )
  
  ################################## Panel B: External transmission (flux shaping) ##################################
  
  # Three scenarios with flux_shape (kappa) values
  scenarios_external <- data.frame(
    name = c("Scenario 1", "Scenario 2", "Scenario 3"),
    flux_shape = c(1.5, 0.001, -1),
    flux_scale = c(1, 1, 1)  # Always 1
  )
  
  # Create combined legend data
  legend_data <- data.frame(
    name = c("Scen. 1", "Scen. 2", "Scen. 3"),
    kappa = scenarios_external$flux_shape
  )
  
  # Create sequence of raw flux values (probability from 0 to 1)
  raw_flux <- seq(0, 1, length.out = 1000)
  
  # Calculate shaped flux for each scenario
  # φ = ν * (1 - e^(-κ*F)) / (1 - e^(-κ))
  results_external <- data.frame()
  
  for (i in 1:nrow(scenarios_external)) {
    kappa <- scenarios_external$flux_shape[i]
    nu <- scenarios_external$flux_scale[i]
    
    # Handle the special case where kappa ≈ 0 (use linear approximation)
    if (abs(kappa) < 1e-6) {
      shaped_flux <- nu * raw_flux  # Linear (no reshaping)
    } else {
      shaped_flux <- nu * (1 - exp(-kappa * raw_flux)) / (1 - exp(-kappa))
    }
    
    results_external <- rbind(results_external, data.frame(
      raw_flux = raw_flux,
      shaped_flux = shaped_flux,
      scenario = scenarios_external$name[i],
      flux_shape = kappa
    ))
  }
  
  results_external$scenario <- factor(results_external$scenario, levels = scenarios_external$name)
  
  # Create Panel B
  panel_B <- ggplot(results_external, aes(x = raw_flux, y = shaped_flux, color = scenario)) +
    geom_line(linewidth = linewidth) +
    
    # Add diagonal reference line (y = x, no reshaping) - make more visible
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "gray30", linewidth = linewidth * 1.2) +
    
    # Scales and labels
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_color_manual(
      values = scenario_colors,
      name = NULL,
      labels = c(
        expression(paste("Scen. 1: ", italic(κ), "=1.5")),
        expression(paste("Scen. 2: ", italic(κ), "=0.0")),
        expression(paste("Scen. 3: ", italic(κ), "=-1.0"))
      )
    ) +
    
    labs(
      x = expression(italic(F)),
      y = expression(italic(φ))
    ) +
    
    # Theme
    theme_classic(base_family = "Georgia") +
    theme(
      axis.title = element_text(size = titlesize, color = 'black'),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.text = element_text(size = textsize, color = 'black'),
      axis.ticks = element_line(color = "black"),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.text = element_text(size = legend_textsize, lineheight = 1.3),
      legend.key.height = unit(0.6, "cm"),
      legend.background = element_rect(fill = "white", color = "gray70", 
                                       linewidth = 0.3),
      legend.margin = margin(5, 5, 5, 5),
      plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
    )
  
  ################################## Combine panels ##################################
  
  combined_plot <- panel_A + panel_B + 
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = titlesize, family = "Georgia"),
          plot.tag.position = c(0.98, 0.85))
  
  # Display plot
  print(combined_plot)
  
  # Save plot
  output_dir <- here("output", "output_figures_tables")
  
  ggsave(
    filename = file.path(output_dir, "transmission_threshold_behavior.png"),
    plot = combined_plot,
    width = fig_width,
    height = fig_height,
    dpi = 300,
    bg = "white"
  )
  
  cat("\n✓ Plot saved to:", file.path(output_dir, "transmission_threshold_behavior.png"), "\n")
