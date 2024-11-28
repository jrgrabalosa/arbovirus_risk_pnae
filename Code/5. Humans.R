library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Població estacional ----------------------------------------------------------
# Data from https://www.idescat.cat/pub/?id=epe&n=9523&geo=com:02

stational_population <- read.csv("Data/Human_data/Població_estacional_Alt_Empordà.csv")

# Rename variables
colnames(stational_population)[colnames(stational_population) == "any"] ="year"
colnames(stational_population)[colnames(stational_population) == "trimestre"] ="trimester"
colnames(stational_population)[colnames(stational_population) == "població_estacional"] ="stational_population"

stational_population <- stational_population %>%
  dplyr::select(year, T1, T2, T3, T4) %>%
  pivot_longer(!c(year),
               names_to = "trimester", values_to = "stational_population")

stational_population$year <- as.factor(stational_population$year)
stational_population$trimester <- as.factor(stational_population$trimester)

stational_population %>%
  ggplot(aes(x = trimester, y = stational_population, fill = year, group = year)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2", name = "Year") +
  labs(x = "Trimester", y = "Stational population", fill = "Year")

# Mean stational population
mean_population <- stational_population %>%
  filter(year %in% c("2015", "2016", "2017", "2018")) %>%
  group_by(trimester) %>%
  summarise(mean_population = mean(stational_population))

mean_population <- mean_population %>%
  mutate(
    initial_week = c(1, 14, 27, 40),
    final_week = c(13, 26, 39, 52)
  )

# Add Aedes albopictus prediction ----------------------------------------------
# Read prediction
pred <- read.csv("Prediction/Aedes_albopictus/pred.csv")

agg_newdata <- pred %>%
  dplyr::select(year, week, pp) %>% 
  group_by(year, week) %>% 
  summarise(pp = mean(pp))

agg_newdata$year <- as.factor(agg_newdata$year)
agg_newdata$week <- as.integer(as.character(agg_newdata$week))

mean_counts <- agg_newdata %>%
  dplyr::select(week, pp) %>% 
  group_by(week) %>%
  summarise(pp = mean(pp))

mean_counts <- mean_counts %>%
  mutate(pp_scaled = pp * 10000)  # Escalar pp para que coincida con el eje primario

# Plot -------------------------------------------------------------------------
albopictus_humans_plot <- ggplot() +
  geom_rect(data = mean_population,
    aes(xmin = initial_week, xmax = final_week, ymin = 0, ymax = mean_population, fill = "Human population"),
    alpha = 0.6, color = "black", size = 0.3) +  # Black border and adjusted opacity for the bars
  geom_line(
    data = mean_counts, 
    aes(x = week, y = pp_scaled, color = "Aedes albopictus"), 
    size = 1.2) +
  scale_y_continuous(
    name = "FTE seasonal human population",
    sec.axis = sec_axis(~ . / 10000, name = "Mean mosquito counts per trap")) +  # Adjusted right axis
  scale_x_continuous(
    breaks = seq(1, max(agg_newdata$week), by = 5)) +  # Define custom X-axis tick marks
  labs(x = "Week", fill = NULL, color = NULL) +
  scale_fill_manual(values = c("Human population" = "grey")) +
  scale_color_manual(values = c("Aedes albopictus" = "#d62728")) +
  theme_minimal(base_size = 14) +  
  theme(
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),  
    axis.title.y.right = element_text(size = 16, margin = margin(l = 10)), 
    axis.text.x = element_text(size = 14, color = "black"),  
    axis.text.y = element_text(size = 14, color = "black"),  
    axis.line = element_line(color = "black", size = 0.5),  
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.89, 0.55),
    legend.text = element_text(size = 14),
    legend.key.size = unit(0.7, "cm"),
    legend.spacing.y = unit(0, "cm"),
    # legend.background = element_rect(fill = "white"),
    legend.box = "vertical") +  # Group items in a single box
  scale_color_manual(
    labels = expression(italic("Aedes albopictus")),
    values = c("Aedes albopictus" = "#d62728")
  )

albopictus_humans_plot

filename <- "Plots/Humans/albopictus_humans.jpg"
ggsave(filename, dpi = 300, width = 15, height = 8, units = "in", type = "jpg", quality = 100)

