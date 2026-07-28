#Model Fit Calculations
library(ggplot2)
library(ggpubr)
library(tidyverse)

#Must read in atx 

###############------------------------------------------------------------------
#Read in all data needed for calculations

                               #River-Wide

#Read in OBSERVED states of River-Wide
obs_river_data <- readRDS(here::here("data/Outputs for Sims and Model Fits/obs_river_data.rds"))
#Extract log-transformed observed data vectors from obs_river_data, 
y_micro <- obs_river_data$microcoleus
y_ana <- obs_river_data$anabaena_cylindrospermum
y_green <- obs_river_data$green_algae
y_nfix <- obs_river_data$other_nfixers
#Put them into a single matrix
y_obs_river <- rbind(y_green, y_micro, y_ana, y_nfix) #Dimensions: Species, time

#Read in LATENT states of River-Wide
allfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_AllVar_predictions.rds"))
bioticfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_Biotic_predictions.rds"))
abioticfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_Abiotic_predictions.rds"))
abioticnonutfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_AbioticNonut_predictions.rds"))

#Read in SIMULATED states of River-Wide
pred.allfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_AllVar.rds"))
pred.bioticfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Biotic.rds"))
pred.abioticfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Abiotic.rds"))
pred.abioticnonutfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_AbioticNoNut.rds"))


                               #Within-Mat

#Read in OBSERVED states of Within-Mat: TM
obs_mat_data <- readRDS(here::here("data/Outputs for Sims and Model Fits/obs_TM_data.rds"))
#Extract log-transformed observed data vectors from obs_mat_data, raw data object
y_ana <- obs_mat_data$Anabaena                  #1
y_epi <- obs_mat_data$`Epithemia Diatoms`       #2
y_geit <- obs_mat_data$Geitlerinema             #3
#Put them into a single matrix
y_obs_mat <- rbind(y_ana, y_epi, y_geit) #Dimensions: Species, time

#Read in LATENT states of Within-Mat: TM
TMfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/WithinMat_Micro_predictions.rds"))

#Read in SIMULATED states of Within-Mat
TM.pred <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/WithinMat_Pred_TM.rds"))

###############------------------------------------------------------------------
#Create lists of the models to iterate through in for loop
river_species <- c("green", "micro", "ana", "nfix")
mat_species <- c("ana", "epi", "geit")

#Create list of all models       <- Later add in output for Toxin models
model_list <- list(
  list(
    model = "allfit",         #The model's name
    category = "River-Wide",  #Whether it used percent cover or microscopy data
    y_obs = y_obs_river,      #Reading in observed values
    post = allfit[["n"]],     #Reading in latent states
    pred = pred.allfit,       #Reading in simulated/predicted values
    species = river_species   #List of species names used in this model
  ),
  list(
    model = "bioticfit",
    category = "River-Wide",
    y_obs = y_obs_river,
    post = bioticfit[["n"]],
    pred = pred.bioticfit,
    species = river_species
  ),
  list(
    model = "abioticfit",
    category = "River-Wide",
    y_obs = y_obs_river,
    post = abioticfit[["n"]],
    pred = pred.abioticfit,
    species = river_species
  ),
  list(
    model = "abioticnonutfit",
    category = "River-Wide",
    y_obs = y_obs_river,
    post = abioticnonutfit[["n"]],
    pred = pred.abioticnonutfit,
    species = river_species
  ),
  list(
    model = "TMfit",
    category = "Within-Mat",
    y_obs = y_obs_mat,
    post = TMfit[["n"]],
    pred = TM.pred,
    species = mat_species
  )
)

###############------------------------------------------------------------------
                            #Calculating model fit indices

#Create an empty dataframe for storing all calculated indices
model.indices <- data.frame(
  model = character(),
  category = character(),
  species = character(),
  metric = character(), #Name of the metric (Bayesian R2, R2, RMSE)
  value = numeric(),    #Value of the metric
  lwr = numeric(),      #Upper confidence interval
  upr = numeric(),      #Lower confidence interval
  stringsAsFactors = FALSE #Do not make categorical variables factors
)

#Set starting counter to 1
counter <- 1

for (m in 1:length(model_list)) {
  
  model <- model_list[[m]]          #Start with the first model of the list
  
  y_obs <- model$y_obs              #Observed data
  posteriors <- exp(model$post)     #For comparing latent values on the same (non-logged) scale as observed and predicted data
  predictives <- model$pred
  species_names <- model$species
  
  iter <- dim(posteriors)[1]     #Number of iterations
  species <- dim(posteriors)[2]  #Number of species/taxa
  time <- dim(posteriors)[3]     #Number of timesteps
  
  #Create metric storage vectors
  BayesR2_vals <- matrix(NA, iter, species) 
  RMSE_vals <- matrix(NA, iter, species) 
  NRMSE_vals <- matrix(NA, iter, species)
  R2_vals <- matrix(NA, iter, species) 
  
  for(s in 1:species) {
    
    #Pull current observation point
    y_obs_s <- y_obs[s, ]
    obs_index <- which(y_obs_s != -99) #Remove weeks where we did not collect field data
    
    for (i in 1:iter) {
      #Bayesian R2      
      yhat <- posteriors[i, s, obs_index] #In the time index position, take out weeks with no field observations
      yobs <- y_obs_s[obs_index]
      
      var_fit <- var(yhat)
      var_res <- var(yobs - yhat)
      
      BayesR2_vals[i, s] <- var_fit / (var_fit + var_res)
      
      #RMSE
      y <- posteriors[i, s, obs_index]
      y_pred <- predictives[i, s, obs_index]
      
      RMSE_vals[i, s] <- sqrt(mean((y - y_pred)^2)) #Calculate RMSE per species iteration
      
      range_y <- max(y) - min(y)
      
      NRMSE_vals[i, s] <- RMSE_vals[i, s] / range_y
      
      #R2      
      R2_vals[i, s] <- cor(y, y_pred)^2 #Calculate R2 per species iteration
      
      #WAIC
      
      
    }
  } 
  #Calculate metrics
  #Consolidate metrics per iter
  metrics <- list(
    "Bayesian R2" = BayesR2_vals,
    "RMSE" = RMSE_vals, 
    "NRMSE" = NRMSE_vals,
    "r2" = R2_vals
  )
  
  #Calculate median and CIs for each metric
  for(j in names(metrics)){
    
    vals <- metrics[[j]]
    
    med <- apply(vals, 2, median)
    ci <- apply(vals, 2, quantile, c(0.025, 0.975))
    
    for(s in 1:species){
      model.indices[counter, ] <- data.frame(
        model = model$model,
        category = model$category,
        species = species_names[s],
        metric = j,
        value = med[s],
        lwr = ci[1,s],
        upr = ci[2,s],
        stringsAsFactors = FALSE
      )
      
      counter <- counter + 1
      
    }
  }
}

###############------------------------------------------------------------------
#Quick visualization of indices so far
str(model.indices)

clean.model.indices <- model.indices %>% 
  dplyr::mutate(species = case_when(species == "ana" ~ "Anabaena",
                                    species == "epi" ~ "Epithemia Diatoms",
                                    species == "geit" ~ "Geitlerinema",
                                    species == "green" ~ "Green Algae",
                                    species == "micro" ~ "Microcoleus",
                                    species == "nfix" ~ "Other N Fixers")) %>% 
  dplyr::mutate(model = case_when(model == "allfit" ~ "All Variables",
                                    model == "bioticfit" ~ "Biotic Only",
                                    model == "abioticfit" ~ "Abiotic Only",
                                    model == "abioticnonutfit" ~ "Abiotic Minus Nutrients",
                                    model == "TMfit" ~ "Target Microcoleus"))

clean.model.indices$model <- factor(  #Manually order model name 
  clean.model.indices$model,
  levels = c("All Variables", "Biotic Only", "Abiotic Only", "Abiotic Minus Nutrients", "Target Microcoleus")
)
clean.model.indices$species <- factor(  #Manually order model name 
  clean.model.indices$species,
  levels = c("Microcoleus", "Anabaena", "Green Algae", "Other N Fixers")
)

#Create a color palette
mycols <- c("brown", "darkolivegreen4", "darkcyan", "darkorange", "#00538A", "#F6926A")
mypal <- palette(mycols)
names(mypal) <- c("Anabaena", "Green Algae", "Microcoleus", 
                  "Other N Fixers", "Epithemia Diatoms", "Geitlerinema")
colScale <- scale_color_manual(values = mypal)
myshap <- c(16, 17, 15, 6, 8, 9)
names(myshap) <- c("Anabaena", "Green Algae", 
                   "Microcoleus", "Other N Fixers", "Epithemia Diatoms", 
                   "Geitlerinema")
shapScale <- scale_shape_manual(values = myshap)

#Plot River-Wide Metrics
ggplot(subset(clean.model.indices, metric %in% c("NRMSE", "RMSE") & category %in% "River-Wide"), aes(x = value, y = species, shape = model, color = model)) +
  facet_wrap(~ metric, scales = "free_x") +
  geom_point(position = position_dodge(width = 0.6), #position_dodge seps species apart
             size = 3) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2,
                 position = position_dodge(width = 0.6)) +
  scale_y_discrete(limits = rev(levels(clean.model.indices$species)[1:4])) + #Reverses order of yaxis
  #colScale + shapScale +
  theme_bw() +
  labs(x = "Metric Value",
       y = "Model Name",
       title = "River-Wide with DIN Goodness-of-Fit",
       shape = "Model",
       color = "Model")


ggplot(subset(subset(clean.model.indices, category %in% "River-Wide"), metric == "RMSE"), aes(x = value, y = model, shape = species, color = species)) +
  geom_point(position = position_dodge(width = 0.6), #position_dodge seps species apart
             size = 3) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2,
                 position = position_dodge(width = 0.6)) +
  scale_y_discrete(limits = rev(levels(clean.model.indices$model)[1:4])) +
  colScale + shapScale +
  coord_cartesian(xlim = c(0, 10)) +
  theme_bw() +
  labs(title = "RMSE",
       x = "Metric value",
       y = "Model",
       shape = "Species",
       color = "Species")

#Plot Within-Mat Metrics
ggplot(subset(clean.model.indices, category %in% "Within-Mat"), aes(x = value, y = model, shape = species, color = species)) +
  facet_wrap(~ metric, scales = "free_x") +
  geom_point(position = position_dodge(width = 0.4), #position_dodge seps species apart
             size = 3) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2,
                 position = position_dodge(width = 0.)) +
  scale_y_discrete(limits = rev(levels(clean.model.indices$model)[5])) + #Reverses order of yaxis
  colScale + shapScale +
  theme_bw() +
  labs(x = "Metric Value",
       y = "Model Name",
       title = "Within-Mat Goodness-of-Fit",
       shape = "Species",
       color = "Species")
