library(readxl)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(lubridate)
library(janitor)
library(vegan)
library(factoextra)
library(tibble)
library(MASS)
library(svglite)
library(pracma)
library(zoo)
library(ggthemes)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# AUC per experiment and subsequent averaging per species ----------------------
# This way we are being more "conservative". If an experiment contains a value 
# above 5, that species will have an AUC > 0 and we will consider it as a reservoir.

# Delete all variables.
rm(list = ls())

# Read data
PNAE_birds <- read_excel("Data/Birds_data/WNV_list/PNAE_birds.xlsx") %>% 
  dplyr::select(latin_name, old_latin_name, family, order, abundance)
viremia_experiments <- read_excel("Data/Birds_data/WNV_list/Review_experimental_infections_birds.xlsx") %>%
  clean_names()

# Editing
viremia_experiments$day <- as.numeric(viremia_experiments$day)
viremia_experiments <- viremia_experiments %>%
  filter(!is.na(day)) %>%
  group_by(curve) %>%
  mutate(
    individuals_tested = max(sample_size_original),
    survival = min(sample_size_survival)/max(sample_size_original)) %>%
  dplyr::select(species, host_family, host_order, country, wnv_strain, virus_genotype, virus_dose, 
                curve, individuals_tested, survival, day, titer) %>%
  arrange(curve, day)

# Add day 0
viremia_experiments <- viremia_experiments %>%
  group_by(curve) %>%
  arrange(curve, day) %>%
  slice(1) %>% # Take the first row of each group (first day)
  mutate(day = 0, titer = 0) %>% # Change the day to 0 and the titer to 0
  bind_rows(viremia_experiments) %>% # Combine with the original dataframe
  arrange(curve, day)

# NAs interpolation
viremia_experiments <- viremia_experiments %>%
  mutate(titer_interpolated = na.approx(titer, day, na.rm = FALSE)) %>%
  filter(!is.na(titer_interpolated))

# Dataframe of viremia experiments filtering PNAE species
viremia_experiments_PNAE <- viremia_experiments %>% 
  filter(species %in% c(PNAE_birds$latin_name)) %>%
  dplyr::select(curve, species, host_family, host_order, day, titer, titer_interpolated)

# Function to calculate AUC for titer values above 5 PFU
calculate_auc <- function(df, threshold = 5) {
  
  # Adjust all titer values: if below threshold, set to threshold
  df_adjusted <- df %>% mutate(titer_adjusted = pmax(titer_interpolated, threshold))
  
  # Calculate the AUC using adjusted titer values
  return(trapz(df_adjusted$day, df_adjusted$titer_adjusted - threshold))
}

# Apply function to each experiment
auc_per_experiment <- viremia_experiments_PNAE %>%
  group_by(species, host_family, host_order, curve) %>%
  summarise(auc = calculate_auc(cur_data()), .groups = 'drop')

# Mean by species, family and order
auc_per_species_1 <- auc_per_experiment %>%
  group_by(species, host_family, host_order) %>%
  summarise(auc_species_1 = mean(auc))

auc_per_family_1 <- auc_per_experiment %>%
  group_by(host_family, host_order) %>%
  summarise(auc_family_1 = mean(auc))

auc_per_order_1 <- auc_per_experiment %>%
  group_by(host_order) %>%
  summarise(auc_order_1 = mean(auc))

# Merge with PNAE list
auc_per_species_1$species_2 <- auc_per_species_1$species
PNAE_birds_WNV <- merge(PNAE_birds, auc_per_species_1, 
                        by.x = c("latin_name", "old_latin_name", "family", "order"),
                        by.y = c("species", "species_2", "host_family", "host_order"), 
                        all.x = TRUE)
PNAE_birds_WNV <- merge(PNAE_birds_WNV, auc_per_family_1, 
                        by.x = c("family", "order"), by.y = c("host_family", "host_order"), 
                        all.x = TRUE)
PNAE_birds_WNV <- merge(PNAE_birds_WNV, auc_per_order_1, 
                        by.x = c("order"), by.y = c("host_order"), 
                        all.x = TRUE)

PNAE_birds_WNV <- PNAE_birds_WNV %>% 
  dplyr::select(latin_name, old_latin_name, family, order, abundance, 
                auc_species_1, auc_family_1, auc_order_1)

rm(auc_per_experiment, auc_per_family_1, auc_per_order_1, auc_per_species_1, 
   PNAE_birds, viremia_experiments)

# Average curve per species and subsequent AUC calculation ---------------------
# First we calculate an average curve, grouping by species and day, 
# and then we calculate the AUC. This way we are more scientific, we have an 
# average curve per species that groups all experiments together. This eliminates outliers.

# Average curves per species/family/order and day + NAs interpolation
curves_species <- viremia_experiments_PNAE %>%
  dplyr::select(curve, species, host_family, host_order, day, titer) %>%
  group_by(species, host_family, host_order, day) %>%
  summarise(titer = if_else(all(is.na(titer)), NA_real_, mean(titer, na.rm = TRUE))) %>%
  mutate(titer_interpolated = na.approx(titer, day, na.rm = FALSE)) %>%
  filter(!is.na(titer_interpolated))

curves_family <- viremia_experiments_PNAE %>%
  dplyr::select(curve, species, host_family, host_order, day, titer) %>%
  group_by(host_family, host_order, day) %>%
  summarise(titer = if_else(all(is.na(titer)), NA_real_, mean(titer, na.rm = TRUE))) %>%
  mutate(titer_interpolated = na.approx(titer, day, na.rm = FALSE)) %>%
  filter(!is.na(titer_interpolated))

curves_order <- viremia_experiments_PNAE %>%
  dplyr::select(curve, species, host_family, host_order, day, titer) %>%
  group_by(host_order, day) %>%
  summarise(titer = if_else(all(is.na(titer)), NA_real_, mean(titer, na.rm = TRUE))) %>%
  mutate(titer_interpolated = na.approx(titer, day, na.rm = FALSE)) %>%
  filter(!is.na(titer_interpolated))

# Apply function to species/family/order curves
auc_per_species_2 <- curves_species %>%
  group_by(species, host_family, host_order) %>%
  summarise(auc_species_2 = calculate_auc(cur_data()), .groups = 'drop') 

auc_per_family_2 <- curves_family %>%
  group_by(host_family, host_order) %>%
  summarise(auc_family_2 = calculate_auc(cur_data()), .groups = 'drop') 

auc_per_order_2 <- curves_order %>%
  group_by(host_order) %>%
  summarise(auc_order_2 = calculate_auc(cur_data()), .groups = 'drop')

# Merge with PNAE list
auc_per_species_2$species_2 <- auc_per_species_2$species
PNAE_birds_WNV <- merge(PNAE_birds_WNV, auc_per_species_2, 
                        by.x = c("latin_name", "old_latin_name", "family", "order"),
                        by.y = c("species", "species_2", "host_family", "host_order"), 
                        all.x = TRUE)
PNAE_birds_WNV <- merge(PNAE_birds_WNV, auc_per_family_2, 
                        by.x = c("family", "order"), by.y = c("host_family", "host_order"), 
                        all.x = TRUE)
PNAE_birds_WNV <- merge(PNAE_birds_WNV, auc_per_order_2, 
                        by.x = c("order"), by.y = c("host_order"), 
                        all.x = TRUE)

PNAE_birds_WNV <- PNAE_birds_WNV %>% 
  dplyr::select(latin_name, old_latin_name, family, order, abundance, 
                auc_species_1, auc_family_1, auc_order_1,
                auc_species_2, auc_family_2, auc_order_2)
writexl::write_xlsx(PNAE_birds_WNV, "Data/Birds_data/WNV_list/PNAE_birds_WNV.xlsx")

# Plots ------------------------------------------------------------------------
# Plot 1: mean curves by family
# Calculate the maximum value of titer_interpolated for each host_family
family_max_titer <- curves_family %>%
  group_by(host_family) %>%
  summarise(max_titer = max(titer_interpolated, na.rm = TRUE)) %>%
  arrange(desc(max_titer))

# Reorder levels of host_family factor based on max_titer values
curves_family <- curves_family %>%
  mutate(host_family = factor(host_family, levels = family_max_titer$host_family))

# Create the plot
mean_curves_family_plot <- curves_family %>%
  filter(!day %in% c("8", "9", "10")) %>%
  ggplot(aes(x = day, y = titer_interpolated, color = host_family)) +
  geom_line() + 
  geom_point(size = 1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", size = 0.5) +  # Horizontal line at threshold (y = 5)
  scale_x_continuous(breaks = seq(0, max(curves_family$day, na.rm = TRUE), by = 1)) +
  labs(title = "Viremia Levels by Family",
       x = "Day postinfection",
       y = "Log PFU/ml serum",
       color = "Family") +
  theme_minimal() +
  theme(legend.position = "right") # Adjust legend position if necessary

mean_curves_family_plot
filename <- "Plots/Birds/mean_curves_family.jpg"
ggsave(filename, dpi = 300, width = 10, height = 6, units = "in", type = "jpg", quality = 100)

rm(curves_family, family_max_titer, mean_curves_family_plot)

# Plot 2: mean curves by order
# Calculate the maximum value of titer_interpolated for each host_order
order_max_titer <- curves_order %>%
  group_by(host_order) %>%
  summarise(max_titer = max(titer_interpolated, na.rm = TRUE)) %>%
  arrange(desc(max_titer))

# Reorder levels of host_order factor based on max_titer values
curves_order <- curves_order %>%
  mutate(host_order = factor(host_order, levels = order_max_titer$host_order))

# Create the plot
mean_curves_order_plot <- curves_order %>%
  filter(!day %in% c("8", "9", "10")) %>%
  ggplot(aes(x = day, y = titer_interpolated, color = host_order)) +
  geom_line() + 
  geom_point(size = 1) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black", size = 0.5) +  # Horizontal line at threshold (y = 5)
  scale_x_continuous(breaks = seq(0, max(curves_order$day, na.rm = TRUE), by = 1)) +
  labs(title = "Viremia Levels by Order",
       x = "Day postinfection",
       y = "Log PFU/ml serum",
       color = "Order") +
  theme_minimal() +
  theme(legend.position = "right") # Adjust legend position if necessary

mean_curves_order_plot
filename <- "Plots/Birds/mean_curves_order.jpg"
ggsave(filename, dpi = 300, width = 10, height = 6, units = "in", type = "jpg", quality = 100)

rm(order_max_titer, mean_curves_order_plot, curves_order)

# Plot 3: experiments + average curves per species
# Preparing dataframe
curves_species$curve <- "mean"
curves_species <- curves_species %>%
  dplyr::select(curve, species, host_family, host_order, day, titer, titer_interpolated)
viremia_experiments_PNAE$curve <- as.character(viremia_experiments_PNAE$curve)
viremia_experiments_PNAE <- rbind(viremia_experiments_PNAE, curves_species)

# Plot
species_curves_plot <- viremia_experiments_PNAE %>%
  ggplot(aes(x = day, y = titer_interpolated, group = curve, color = ifelse(curve == "mean", "mean", "experiment"))) +
  geom_line() +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 10, by = 1), limits = c(0, 10)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10)) +
  scale_color_manual(values = c("mean" = "black", "experiment" = "darkolivegreen3")) +
  facet_wrap(~species, scales = "free") +
  labs(
    # title = "Experiments and Average Viremia Curves per Species",
    x = "Day postinfection",
    y = "Log PFU/ml serum") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face = "italic", size = 10),
        axis.text = element_text(size = 10),  # Tamaño del texto en los ejes
        axis.title = element_text(size = 11))  # Tamaño de los títulos de los ejes

species_curves_plot
filename <- "Plots/Birds/species_curves.jpg"
ggsave(filename, dpi = 300, width = 10, height = 6, units = "in", type = "jpg", quality = 100)

rm(curves_species, species_curves_plot, viremia_experiments_PNAE, filename, calculate_auc)

# # Per calcular el nombre d'experiments i l'interval de confiança
# viremia_experiments_PNAE <- viremia_experiments_PNAE %>%
#   group_by(species) %>%
#   mutate(n_experiments = n_distinct(curve)) %>%
#   ungroup() 
# 
# viremia_summary <- viremia_experiments_PNAE %>%
#   group_by(species) %>%
#   mutate(
#     n_experiments = n_distinct(curve) 
#   ) %>%
#   ungroup() %>%
#   group_by(species, day, n_experiments) %>%
#   summarise(
#     mean_titer = mean(titer_interpolated, na.rm = TRUE),
#     sd_titer = sd(titer_interpolated, na.rm = TRUE)
#   )

################################################################################
# Titer-infectiousness curve
################################################################################

# Data from Komar, N., Langevin, S., Hinten, et al. (2003). Experimental infection 
# of North American birds with the New York 1999 strain of West Nile virus. 
# Emerging infectious diseases, 9(3), 311–322. https://doi.org/10.3201/eid0903.020628
# "To produce these data, we used a threshold level of infectious viremia 
# of 105.0 PFU/mL serum and estimated infectiousness of each bird’s viremia levels 
# from a standard curve for infection of C. pipiens as a function of viremic titer 
# derived from Turell et al. (13) (Appendix D)."

# Delete all variables.
rm(list = ls())

# Read data
titer_infectiousness_curve <- read_excel("Data/Birds_data/WNV_list/titer_infectiousness_curve.xlsx")
# extra_data <- data.frame(Viremia = seq(0, 4.9, by = 0.1), Infectiousness = 0)
# titer_infectiousness_curve <- rbind(titer_infectiousness_curve, extra_data) %>%
#   arrange(Viremia)

titer_infectiousness_plot <- ggplot(data = titer_infectiousness_curve, aes(x = Viremia, y = Infectiousness)) +
  geom_line() +
  labs(x = "Viremia (log PFU/ml serum)", y = "Infectiousness",
       title = "Standard curve relating viremia titer to infectiousness for Culex pipiens") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 12, by = 1)) +  # Etiquetas en ejes x como enteros
  scale_y_continuous(breaks = seq(0, 0.7, by = 0.1))  # Etiquetas en ejes y en decimales

titer_infectiousness_plot
# filename <- "Plots/Birds/titer_infectiousness_plot.jpg"
# ggsave(filename, dpi = 300, width = 10, height = 6, units = "in", type = "jpg", quality = 100)

