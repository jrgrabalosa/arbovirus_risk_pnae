# Arboviral disease risk in a Mediterranean wetland area

Here, we assess the potential risk of arboviral transmission in a high-risk, non-endemic Mediterranean region of northern Catalonia, analyzing mosquito vector populations (Aedes albopictus, Culex pipiens, Culex modestus and Culex theileri) and West Nile virus (WNV) avian hosts.

# Installation

Make sure to have R 4.3.2 (or later) and, optionally, RStudio installed for an enhanced development environment.

1. Clone the repository:
git clone https://github.com/username/project.git
cd project

2. Install required packages: Open R or RStudio and run the following code to install the required packages:

packages <- c("readxl", "writexl", "tidyverse", "lubridate", "dplyr", "ggplot2", "purrr", "zoo", "parallelly", "parallel", "janitor", "pollen", "lme4", "DHARMa", "glmmTMB", "performance", "car", "MuMIn", "corrplot", "cmdstanr", "brms", "rstanarm", "loo", "tidybayes", "RColorBrewer", "sf", "raster", "viridis", "ggrepel", "vegan", "factoextra", "tibble", "MASS", "svglite", "pracma", "ggthemes")

packages_to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(packages_to_install)) install.packages(packages_to_install) #Install packages that are not already installed

lapply(packages, library, character.only = TRUE) #Load packages

3. Run the project: Open the main script file (e.g., main.R) in R or RStudio and execute it:
source("main.R")

# Project Directory Structure

The project directory PNAE is organized as follows:

```plaintext
PNAE/
├── Code/                  		# Contains all R scripts and code for analysis
│   │
│   ├── 1_Aedes_albopictus.R    	# Dataframe preparation, models, and predictions for Aedes albopictus data
│   ├── 1_Culex.R               	# Dataframe preparation, models, and predictions for Culex spp. data
│   ├── 2_Birds_experiments_review.R   	# Host competence calculations using experimental infection data
│   ├── 3_Birds.R               	# WNV host spatial analysis
│   ├── 4_Risk.R                	# Overlapping mosquitoes and birds
│   ├── 5_Humans.R              	# Human population dynamics
│   ├── 6_Cross_validation.R    	# Model cross-validation
│   │	
├── Data/                 		# Contains data files and output datasets
│   │ 
│   ├── Birds_data/       		# Data from the Atlas of Nesting Birds of Catalonia and review of experimental WNV infection studies
│   │   ├── Atlas_data/
│   │   ├── Outputs/
│   │   └── WNV_list/
│   ├── Human_data/       		# Data related to human population and demographics
│   ├── Land_covers/      		# Land cover and environmental data
│   │   └── UsosCobertes_2017/
│   ├── Meteo_data/       		# Meteorological data
│   ├── Mosquito_data/    		# SCM mosquito surveillance and output databases
│   │   ├── Aedes_albopictus/
│   │   ├── Culex/
│   │   └── Mosquito_Surveillance_SCM/
│   └── Risk_data/        		# Risk assessment data
│
├── Models/               		# Statistical models used in the project
│   │
│   ├── Aedes_albopictus/
│   │   ├── brms_models/  		# Bayesian models for Aedes albopictus
│   │   └── glm_models/   		# GLM models for Aedes albopictus
│   └── Culex/
│       ├── brms_models/  		# Bayesian models for Culex species
│       └── glm_models/   		# GLM models for Culex species
│
├── Plots/                		# Plots and visualizations generated from data
│   │
│   ├── Aedes_albopictus/
│   ├── Birds/
│   └── Culex/
│
└── Predictions/          		# Contains space-time prediction files
    │
    ├── Aedes_albopictus/
    └── Culex/



