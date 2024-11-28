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

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Working with the Atlas data --------------------------------------------------
# Read data
atlas_data <- read_delim("Data/Birds_data/Atlas_data/Estimes_WNV.csv", delim = ";")
colnames(atlas_data)[colnames(atlas_data) == "latin_name"] <- "atlas_name"
area_interes <- read_excel("Data/Birds_data/area_interes.xlsx")
PNAE_birds_WNV <- read_excel("Data/Birds_data/WNV_list/PNAE_birds_WNV.xlsx")
correspondencia <- PNAE_birds_WNV %>% dplyr::select(latin_name, old_latin_name)

# Filter by study site (PNAE)
atlas_data <- merge(atlas_data, area_interes, by.x = "utm_1x1", by.y = "UTM1X1") %>%
  clean_names() %>%
  dplyr::select(atlas_name, utm_1x1, zona, minim, maxim) # filter by area

# Fix names
atlas_names <- unique(atlas_data$atlas_name)
atlas_names_df <- as.data.frame(atlas_names)
latin_names <- unique(PNAE_birds_WNV$latin_name)
old_latin_names <- unique(PNAE_birds_WNV$old_latin_name)

result_list <- list()

for(i in 1:nrow(atlas_names_df)) {
  atlas_name <- atlas_names_df$atlas_name[i]
  
  # Search for a match in 'old_latin_name'
  match_old <- correspondencia[correspondencia$old_latin_name == atlas_name, ]
  
  # Search for a match in 'latin_name'
  match_new <- correspondencia[correspondencia$latin_name == atlas_name, ]
  
  # If there is a match in 'old_latin_name'
  if(nrow(match_old) > 0) {
    result_list[[i]] <- data.frame(
      atlas_name = atlas_name,
      old_latin_name = match_old$old_latin_name,
      latin_name = match_old$latin_name
    )
  }
  # If there is a match in 'latin_name'
  else if(nrow(match_new) > 0) {
    result_list[[i]] <- data.frame(
      atlas_name = atlas_name,
      old_latin_name = match_new$old_latin_name,
      latin_name = match_new$latin_name
    )
  }
}

correspondencia <- do.call(rbind, result_list)
rm(atlas_names_df, match_new, match_old, result_list, atlas_name, atlas_names,
   latin_names, old_latin_names, i)

PNAE_birds_WNV_atlas <- merge(PNAE_birds_WNV, correspondencia, by = c("latin_name", "old_latin_name")) %>%
  dplyr::select(atlas_name, family, order, abundance, auc_species_2)

# Add viremia
atlas_data <- merge(atlas_data, PNAE_birds_WNV_atlas, by = c("atlas_name"), all.x = TRUE) %>%
  dplyr::select(atlas_name, family, order, utm_1x1, zona, minim, maxim, auc_species_2)

# Add weight value (used for further calculations)
atlas_data$weight <- 1

# Host capacity per species
host_capacity <- atlas_data %>%
  filter(!is.na(auc_species_2)) %>%
  group_by(atlas_name, auc_species_2) %>%
  summarise(abundance = sum(maxim)) %>%
  group_by(atlas_name, auc_species_2, abundance) %>%
  summarise(host_capacity = auc_species_2*abundance) %>%
  arrange(desc(host_capacity),desc(abundance))

rm(host_capacity)

# Reservoirs and non-reservoirs abundance
atlas_data_reservoir_species <- atlas_data %>% filter(auc_species_2 > 0)
atlas_data_non_reservoir_species <- atlas_data %>% filter(auc_species_2 == 0)

# Cell calculations ------------------------------------------------------------
cell_pv <- atlas_data %>% 
  group_by(utm_1x1) %>% 
  summarise(number_of_species = sum(weight),
            predicted_abundance = sum(maxim))

# Cell calculations (reservoir species) 
cell_pv_rs <- atlas_data_reservoir_species %>%
  group_by(utm_1x1) %>%
  summarise(number_of_species_rs = sum(weight),
            predicted_abundance_rs = sum(maxim),
            host_capacity_rs = sum(maxim*auc_species_2))

cell_pv <- merge(cell_pv, cell_pv_rs, by = "utm_1x1")

# Cell calculations (non-reservoir species)
cell_pv_nrs <- atlas_data_non_reservoir_species %>%
  group_by(utm_1x1) %>%
  summarise(number_of_species_nrs = sum(weight),
            predicted_abundance_nrs = sum(maxim))

cell_pv <- merge(cell_pv, cell_pv_nrs, by = "utm_1x1")

cell_pv$rs_ratio <- (cell_pv$predicted_abundance_rs)/
  (cell_pv$predicted_abundance_nrs)

# write.csv(cell_pv, "Data/Birds_data/Outputs/birds_UTM1x1.csv")
rm(atlas_data_non_reservoir_species, atlas_data_reservoir_species, 
   cell_pv_nrs, cell_pv_rs)

# Plots
ggplot(cell_pv, aes(x = number_of_species , y = rs_ratio)) +
  geom_point() + 
  geom_smooth(method = "loess", se = FALSE) + 
  # geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) + 
  labs(x = "Number of Species", 
       y = "Reservoirs/Non-reservoirs Ratio", 
       title = "Relationship Between Reservoirs/Non-reservoirs Ratio and Biodiversity") + 
  theme_minimal() 

ggplot(cell_pv, aes(x = predicted_abundance , y = rs_ratio)) +
  geom_point() + 
  geom_smooth(method = "loess", se = FALSE) + 
  # geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) + 
  labs(x = "Birds Abundance", 
       y = "Reservoirs/Non-reservoir Abundance", 
       title = "Relationship Between Reservoirs/Non-reservoirs Ratio and Birds Abundance") + 
  theme_minimal()

# Land uses UTM1x1 -------------------------------------------------------------
land_covers <- read_csv("Data/Land_covers/landcovers_utm1x1.csv") %>%
  group_by(UTM1X1, DN) %>%
  summarise(area = sum(area)) %>%
  pivot_wider(names_from = DN, values_from = area, values_fill = 0)

land_covers$area_total <- rowSums(land_covers[, -1])
land_covers <- land_covers %>%
  mutate(across(1:20, ~ . / area_total)) %>%
  dplyr::select(-area_total)

land_covers$urban_area <- land_covers$`4` + land_covers$`5` + land_covers$`6`
land_covers$wetlands <- land_covers$`10` 
land_covers$forests <- land_covers$`15` + land_covers$`16` + land_covers$`17`
land_covers$dry_crops <- land_covers$`19`
land_covers$irrigated_crops <- land_covers$`20`
land_covers$herbaceous_crops <- land_covers$`19` + land_covers$`20`
land_covers$rice_fields <- land_covers$`21`

land_covers <- land_covers %>%
  dplyr::select(UTM1X1, urban_area, wetlands, forests, dry_crops,
                irrigated_crops, herbaceous_crops, rice_fields)

# RDA --------------------------------------------------------------------------
# Preparing the RDA dataframe (species and land covers in columns)
rda_df <- atlas_data %>%
  dplyr::select(utm_1x1, atlas_name, maxim) %>%
  pivot_wider(names_from = atlas_name, values_from = maxim, values_fill = 0)
rda_df <- merge(rda_df, land_covers, by.x = "utm_1x1", by.y = "UTM1X1")

# Species matrix
species <- as.matrix(rda_df[,2:93]) #matrix

# RDA model
rda_model <- rda(species ~ urban_area + wetlands + irrigated_crops + dry_crops + rice_fields, 
                 data = rda_df, 
                 scale = TRUE, 
                 na.rm = TRUE) #RDA model, land uses as explanatory variables

# Custom RDA plot
# Extract % explained by the first 2 axes
perc <- round(100*(summary(rda_model)$cont$importance[2, 1:2]), 2)

# Extract scores - these are coordinates in the RDA space
sc_sp <- scores(rda_model, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(rda_model, display="bp", choices=c(1, 2), scaling=1)

# Convert to data frames
sc_sp_df <- as.data.frame(sc_sp)
sc_bp_df <- as.data.frame(sc_bp)
sc_sp_df$label <- rownames(sc_sp)
sc_bp_df$label <- c("Urban area", "Wetlands", "Irrigated crops", "Dry crops", "Rice fields")

# Scale arrow coordinates
arrow_scale <- 6.5

# Plot
rda_plot <- ggplot() +
  geom_segment(data = sc_bp_df, aes(x = 0, y = 0, xend = RDA1 * arrow_scale, yend = RDA2 * arrow_scale), 
               arrow = arrow(length = unit(0.3, "cm")), color = "red3", size = 0.6) +
  geom_point(data = sc_sp_df, aes(x = RDA1, y = RDA2), 
             shape = 22, color = "black", fill = "#f2bd33", size = 2, alpha = 0.75) +
  geom_text_repel(data = sc_bp_df, aes(x = RDA1 * arrow_scale, y = RDA2 * arrow_scale, label = label), 
                  color = "red3", size = 3, fontface = "bold") +
  geom_text_repel(data = sc_sp_df, aes(x = RDA1, y = RDA2, label = label), 
                  color = "grey40", size = 2.5, fontface = "bold") +
  xlim(-3, 3) + ylim(-3, 3) +
  labs(title = "Habitat Preferences in Bird Species", 
       x = paste0("RDA1 (", perc[1], "%)"), y = paste0("RDA2 (", perc[2], "%)")) +
  theme_bw()

print(rda_plot)
ggsave(filename = "Plots/Birds/rda_plot.svg", plot = rda_plot, dpi = 300, width = 15, height = 8, units = "in")

# Species coordinates
rda_plot <- plot(rda_model)
rda_axis <- as.data.frame(rda_plot$species) #dataframe of the species and its values on the RDA axes

rm(rda_model, rda_plot, sc_bp, sc_bp_df, sc_sp, sc_sp_df, species, arrow_scale, perc)

# KMEANS -----------------------------------------------------------------------
fviz_nbclust(rda_axis, kmeans, method = "wss") #elbow method
set.seed(123)
k <- 3  # Number of clusters
kmeans_rda <- kmeans(rda_axis, centers = k)
kmeans_plot <- fviz_cluster(kmeans_rda, rda_axis, 
                            ellipse.type = "norm", # Type of ellipse
                            stand = TRUE, # Standardize variables
                            axes = c(1, 2), # Define axes
                            repel = TRUE, # Avoid overlapping
                            show.clust.cent = FALSE, # Don't show cluster centers
                            shape = 20, # Point shape: solid circle
                            pointsize = 1, # Point size
                            labelsize = 8, # Label size
                            main = "Bird species clustering based on habitat selection", # Title
                            legend.title = "Clusters", # Legend title
                            font.main = c(14, "bold"), # Title size
                            font.x = 14, # Axis size
                            font.y = 14, # Axis size
                            font.legend = 14, # Legend size
                            palette = "Dark2")+ # Palette 
  guides(fill = guide_legend(override.aes = list(label = ""))) + # Don't show legend labels
  theme(legend.text = element_text(size = 14), # Legend text size
        legend.key.size = unit(2.5, "lines")) + # Legend key size
  theme_bw()

print(kmeans_plot)
ggsave(filename = "Plots/Birds/kmeans_plot.svg", plot = kmeans_plot, dpi = 300, width = 15, height = 8, units = "in")

# Kmeans analysis
rda_axis$cluster <- kmeans_rda$cluster

for (i in 1:length(rda_axis$cluster)) {
  if (rda_axis$cluster[i] == 1) rda_axis$cluster_name[i] = "Herbaceous crops"
  if (rda_axis$cluster[i] == 2) rda_axis$cluster_name[i] = "Urban area"
  if (rda_axis$cluster[i] == 3) rda_axis$cluster_name[i] = "Wetlands and rice fields"}

rda_axis <- rownames_to_column(rda_axis, var = "latin_name")

rm(kmeans_plot, kmeans_rda, modelo, i, k)
