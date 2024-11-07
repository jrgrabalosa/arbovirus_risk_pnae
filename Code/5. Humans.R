library(tidyverse)
library(ggplot2)

# Delete all variables.
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Tourist data -----------------------------------------------------------------
municipis <- read.csv("Data/Human_data/mun.csv")
tourists <- read.csv("Data/Human_data/slc_tourist.csv")

# Dataframe with tourism data
tourists_mun <- merge(municipis,tourists,by="mun")
tourists_mun$year <- as.character(tourists_mun$year)
tourists_mun <- tourists_mun %>% filter(com == "2") # filter by region "Alt Empordà"

tourists_mun <- tourists_mun %>%
  group_by(mun, name, longitude, latitude, year) %>%
  summarise(acc = sum(f_tourist_acc),
            beds = sum(f_tourist_beds)) # calculation of acc and beds for each municipality and year

tourists_comarca <- tourists_mun %>%
  group_by(year) %>%
  summarise(acc = sum(acc),
            beds = sum(beds)) # calculation of acc and beds for the whole region

# Població estacional ----------------------------------------------------------
# Data from https://www.idescat.cat/pub/?id=epe&n=9523&geo=com:02

població_ETCA <- read.csv("Data/Human_data/Població_estacional_Alt_Empordà.csv")
població_estacional_ETCA <- població_ETCA %>%
  dplyr::select(any, T1, T2, T3, T4) %>%
  pivot_longer(!c(any),
               names_to = "trimestre", values_to = "població_estacional")

població_estacional_ETCA$any <- as.factor(població_estacional_ETCA$any)
població_estacional_ETCA$trimestre <- as.factor(població_estacional_ETCA$trimestre)

plot <- població_estacional_ETCA %>%
  ggplot(aes(x = trimestre, y = població_estacional, fill = any, group = any)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2", name = "Year") +
  labs(x = "Trimestre", y = "Població estacional ETCA", fill = "Year")

plot


