library(readxl)
library(brms)
library(tidyverse)
library(cmdstanr) # Solo si estás usando cmdstanr como backend
library(bayesplot) # Opcional para visualización

# Delete all variables
rm(list = ls())

# Set working directory
setwd("/home/julia/Documentos/PhD/Data_analysis/PNAE")

# Read excel
culex_df <- read_excel("Data/Mosquito_data/Culex/culex_df.xlsx")

# Convert variables
culex_df$trap_name <- as.factor(culex_df$trap_name)
culex_df$AMT <- as.factor(culex_df$AMT)
culex_df$year <- factor(culex_df$year)
culex_df$month <- factor(culex_df$month, ordered = TRUE, levels = as.character(1:12))
culex_df$week <- factor(culex_df$week, ordered = TRUE, levels = as.character(1:53))
culex_df$trap_type <- as.factor(culex_df$trap_type)
culex_df$land_use <- as.factor(culex_df$land_use)

# Función para calcular la cross-validation con el modelo brms
cros_val <- function(culex_df, formula, family, nam = "variables", n = 10) {
  
  # Dataframe vacío para almacenar resultados
  cv_all <- data.frame()
  
  for (i in 1:n) {
    
    cat("Iteración:", i, "\n")
    
    # Dividir datos en conjuntos de entrenamiento y prueba para cada iteración
    sample <- sample(c(TRUE, FALSE), nrow(culex_df), replace = TRUE, prob = c(0.8, 0.2))
    train <- culex_df[sample, ]
    test <- culex_df[!sample, ]
    
    # Ajuste del modelo usando brm y la fórmula especificada
    model <- brm(
      formula = formula,
      data = train,
      family = family,
      prior = set_prior("cauchy(0,2.5)", class = "b"),
      iter = 7000,
      chains = 4,
      cores = 4,
      backend = "cmdstanr",
      control = list(adapt_delta = 0.999),
      silent = TRUE
    )
    
    # Predicción en el conjunto de prueba
    pp <- round(predict(model, newdata = test, re_formula = NA)[, "Estimate"], 0)
    
    # Cálculo de métricas de evaluación
    cv <- data.frame(
      model = nam,
      AIC = AIC(model),
      R2 = cor(pp, test$females)^2, # R2 calculado manualmente como el cuadrado de la correlación entre observados y predichos
      RMSE = sqrt(mean((pp - test$females)^2)), # Raíz del error cuadrático medio
      MAE = mean(abs(pp - test$females)), # Error absoluto medio
      pred.error = sqrt(mean((pp - test$females)^2)) / mean(test$females) # Error de predicción
    )
    
    # Agregar resultados a cv_all
    cv_all <- rbind(cv_all, cv)
    
  }
  
  return(cv_all)
}

# Llamada a la función con tus datos y modelo
cv_all <- cros_val(
  culex_df = culex_df,
  formula = bf(females ~ scale(TM) + scale(HRM7) + scale(PPT21) +
                 (1|year:trap_name) + (1|land_use)),
  family = zero_inflated_negbinomial(link = "log"),
  nam = "Culex ZINB Model",
  n = 10
)
