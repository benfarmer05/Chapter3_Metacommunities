# .rs.restartR(clean = TRUE)
rm(list=ls())

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

linewidths = 0.75
symbsizes = 1.3
titlesize = 10
textsize = 9

############################ SIR MODEL FUNCTION ############################

SIR.multi = function(t, y, p) {
  S.LS = y[1]; I.LS = y[2]; R.LS = y[3]
  S.MS = y[4]; I.MS = y[5]; R.MS = y[6]
  S.HS = y[7]; I.HS = y[8]; R.HS = y[9]
  
  with(as.list(p), {
    P = (I.LS + I.MS + I.HS)
    
    dS.LS.dt = -b.LS * S.LS * P / N.LS
    dI.LS.dt = b.LS * S.LS * P / N.LS - g.LS * I.LS
    dR.LS.dt = g.LS * I.LS
    
    dS.MS.dt = -b.MS * S.MS * P / N.MS
    dI.MS.dt = b.MS * S.MS * P / N.MS - g.MS * I.MS
    dR.MS.dt = g.MS * I.MS
    
    dS.HS.dt = -b.HS * S.HS * P / N.HS
    dI.HS.dt = b.HS * S.HS * P / N.HS - g.HS * I.HS
    dR.HS.dt = g.HS * I.HS
    
    return(list(c(dS.LS.dt, dI.LS.dt, dR.LS.dt, dS.MS.dt, dI.MS.dt, dR.MS.dt, 
                  dS.HS.dt, dI.HS.dt, dR.HS.dt), P = P))
  })
}

##################### FUNCTION TO RUN ONE SCENARIO ####################

run_one_scenario = function(cover_LS, cover_MS, cover_HS, use_cover_scale = TRUE) {
  
  # Prevent division by zero
  epsilon = 1e-6
  cover_LS = max(cover_LS, epsilon)
  cover_MS = max(cover_MS, epsilon)
  cover_HS = max(cover_HS, epsilon)
  
  # Parameters for THIS scenario only
  beta_LS = 0.03
  beta_MS = 0.14
  beta_HS = 2.08
  gamma_LS = 0.05
  gamma_MS = 0.55
  gamma_HS = 3.33
  
  # Scale to 0-1 if requested
  if (use_cover_scale) {
    N_LS = cover_LS / 100
    N_MS = cover_MS / 100
    N_HS = cover_HS / 100
    initial_infected = 0.01 / 100
  } else {
    N_LS = cover_LS
    N_MS = cover_MS
    N_HS = cover_HS
    initial_infected = 0.01
  }
  
  # Seed infection: HS -> MS -> LS priority
  if (cover_HS > epsilon * 10) {
    I0_LS = 0; I0_MS = 0; I0_HS = initial_infected
  } else if (cover_MS > epsilon * 10) {
    I0_LS = 0; I0_MS = initial_infected; I0_HS = 0
  } else {
    I0_LS = initial_infected; I0_MS = 0; I0_HS = 0
  }
  
  # Initial state
  state0 = c(S.LS = N_LS - I0_LS, I.LS = I0_LS, R.LS = 0,
             S.MS = N_MS - I0_MS, I.MS = I0_MS, R.MS = 0,
             S.HS = N_HS - I0_HS, I.HS = I0_HS, R.HS = 0)
  
  # Parameters for THIS scenario
  parms = list(b.LS = beta_LS, b.MS = beta_MS, b.HS = beta_HS,
               g.LS = gamma_LS, g.MS = gamma_MS, g.HS = gamma_HS,
               N.LS = N_LS, N.MS = N_MS, N.HS = N_HS)
  
  # Run model
  times = seq(0, 400, by = 0.1)
  result = as.data.frame(ode(y = state0, times = times, func = SIR.multi, parms = parms))
  
  # Add metadata
  total = cover_LS + cover_MS + cover_HS
  result$scenario = sprintf("L:%.1f M:%.1f H:%.1f (%.1f%%)", 
                            cover_LS, cover_MS, cover_HS, total)
  result$scale_type = ifelse(use_cover_scale, "Cover", "Raw")
  
  return(result)
}

############################# RUN ALL SCENARIOS ##############################

scenario1 = run_one_scenario(3.64, 15.7, 4.97, TRUE)
scenario1$scenario = "Nearshore (cover)"

scenario2 = run_one_scenario(3.64, 15.7, 4.97, FALSE)
scenario2$scenario = "Nearshore (raw)"

scenario3 = run_one_scenario(10, 5, 2, TRUE)
scenario4 = run_one_scenario(5, 2.5, 1, TRUE)
scenario5 = run_one_scenario(2.5, 1.25, 3, TRUE)
scenario6 = run_one_scenario(0.5, 5, 1, TRUE)
scenario7 = run_one_scenario(0.25, 2.5, 0.5, TRUE)
scenario8 = run_one_scenario(0.125, 1.25, 0.25, TRUE)
scenario9 = run_one_scenario(0.1, 1, 0.02, TRUE)

all_results = bind_rows(scenario1, scenario2, scenario3, scenario4, scenario5, 
                        scenario6, scenario7, scenario8, scenario9)

# ################# PLOT 1: Infected Only ##################
# 
# plot_data_infected = all_results %>%
#   select(time, I.LS, I.MS, I.HS, scenario) %>%
#   pivot_longer(cols = c(I.LS, I.MS, I.HS), names_to = "Group", values_to = "Infected") %>%
#   mutate(Group = factor(Group, levels = c("I.LS", "I.MS", "I.HS"),
#                         labels = c("Low", "Moderate", "High")))
# 
# p1 = ggplot(plot_data_infected, aes(x = time, y = Infected, color = Group)) +
#   geom_line(linewidth = linewidths) +
#   facet_wrap(~ scenario, scales = "free_y", ncol = 3) +
#   scale_color_manual(values = c("Low" = "#1E90FF", "Moderate" = "#FFD700", "High" = "#FF1493")) +
#   labs(x = "Day of outbreak", y = "Percent cover", color = "") +
#   theme_classic(base_family = "Georgia") +
#   theme(legend.position = "bottom", axis.title = element_text(size = titlesize),
#         axis.text = element_text(size = textsize), strip.text = element_text(size = 8),
#         legend.text = element_text(size = textsize), legend.title = element_blank(),
#         legend.key.height = unit(0, "cm"), legend.key.size = unit(0.3, "cm"))
# 
# print(p1)
# 
################# PLOT 2: Full SIR ##################

plot_data_full = all_results %>%
  select(time, S.LS, I.LS, R.LS, S.MS, I.MS, R.MS, S.HS, I.HS, R.HS, scenario) %>%
  pivot_longer(cols = -c(time, scenario), names_to = "Compartment_Group", values_to = "Value") %>%
  separate(Compartment_Group, into = c("Compartment", "Group"), sep = "\\.") %>%
  mutate(Group = factor(Group, levels = c("LS", "MS", "HS"), 
                        labels = c("Low", "Moderate", "High")),
         Compartment = factor(Compartment, levels = c("S", "I", "R"),
                              labels = c("Susceptible", "Infected", "Removed")))

p2 = ggplot(plot_data_full, aes(x = time, y = Value, color = Group, shape = Compartment)) +
  geom_line(linewidth = linewidths) +
  facet_wrap(~ scenario, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Low" = "#1E90FF", "Moderate" = "#FFD700", "High" = "#FF1493")) +
  scale_shape_manual(values = c("Susceptible" = 16, "Infected" = 17, "Removed" = 15), name = "") +
  labs(x = "Day of outbreak", y = "Percent cover", color = "") +
  theme_classic(base_family = "Georgia") +
  theme(legend.position = "bottom", axis.title = element_text(size = titlesize),
        axis.text = element_text(size = textsize), strip.text = element_text(size = 8),
        legend.text = element_text(size = textsize), legend.title = element_blank(),
        legend.key.height = unit(0, "cm"), legend.key.size = unit(0.3, "cm"))

print(p2)

# ################# PLOT 3: Cover vs Raw Comparison ##################
# 
# comparison_data = bind_rows(scenario1, scenario2) %>%
#   select(time, S.LS, I.LS, R.LS, S.MS, I.MS, R.MS, S.HS, I.HS, R.HS, scale_type) %>%
#   pivot_longer(cols = -c(time, scale_type), names_to = "Compartment_Group", values_to = "Value") %>%
#   separate(Compartment_Group, into = c("Compartment", "Group"), sep = "\\.") %>%
#   mutate(Group = factor(Group, levels = c("LS", "MS", "HS"), 
#                         labels = c("Low", "Moderate", "High")),
#          Compartment = factor(Compartment, levels = c("S", "I", "R"),
#                               labels = c("Susceptible", "Infected", "Removed")))
# 
# p3 = ggplot(comparison_data, aes(x = time, y = Value, color = Group, shape = Compartment)) +
#   geom_line(linewidth = linewidths) +
#   facet_wrap(~ scale_type, scales = "free_y", ncol = 2) +
#   scale_color_manual(values = c("Low" = "#1E90FF", "Moderate" = "#FFD700", "High" = "#FF1493")) +
#   scale_shape_manual(values = c("Susceptible" = 16, "Infected" = 17, "Removed" = 15), name = "") +
#   labs(x = "Day of outbreak", y = "Tissue amount", color = "") +
#   theme_classic(base_family = "Georgia") +
#   theme(legend.position = "bottom", axis.title = element_text(size = titlesize),
#         axis.text = element_text(size = textsize), strip.text = element_text(size = 8),
#         legend.text = element_text(size = textsize), legend.title = element_blank(),
#         legend.key.height = unit(0, "cm"), legend.key.size = unit(0.3, "cm"))
# 
# print(p3)
# 
# ########################## SUMMARY ############################
# 
# cat("\n=== Model Summary ===\n")
# cat("✓ Each scenario is completely independent\n")
# cat("✓ Infection seeds in HS -> MS -> LS priority\n")
# cat("✓ Zero cover handled with epsilon (1e-6)\n")
# cat("✓ Cover scale (0-1) vs raw abundance comparison included\n")