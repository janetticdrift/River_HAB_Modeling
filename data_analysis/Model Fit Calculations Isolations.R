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
nitratefit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_Nitrate_predictions.rds"))
phosfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_Phos_predictions.rds"))
ammoniumfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_Ammonium_predictions.rds"))
DINfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_DIN_predictions.rds"))

#Read in SIMULATED states of River-Wide
pred.nitratefit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Nitrate.rds"))
pred.phosfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Phosphate.rds"))
pred.ammoniumfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Ammonium.rds"))
pred.DINfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_DIN.rds"))

###############------------------------------------------------------------------
#Create lists of the models to iterate through in for loop
river_species <- c("green", "micro", "ana", "nfix")

#Create list of all models       <- Later add in output for Toxin models
model_list <- list(
  list(
    model = "nitratefit",         #The model's name
    category = "River-Wide",  #Whether it used percent cover or microscopy data
    y_obs = y_obs_river,      #Reading in observed values
    post = nitratefit[["n"]],     #Reading in latent states
    pred = pred.nitratefit,       #Reading in simulated/predicted values
    species = river_species   #List of species names used in this model
  ),
  list(
    model = "phosfit",
    category = "River-Wide",
    y_obs = y_obs_river,
    post = phosfit[["n"]],
    pred = pred.phosfit,
    species = river_species
  ),
  list(
    model = "ammoniumfit",
    category = "River-Wide",
    y_obs = y_obs_river,
    post = ammoniumfit[["n"]],
    pred = pred.ammoniumfit,
    species = river_species
  ),
  list(
    model = "DINfit",
    category = "River-Wide",
    y_obs = y_obs_river,
    post = DINfit[["n"]],
    pred = pred.DINfit,
    species = river_species
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
                                    species == "green" ~ "Green Algae",
                                    species == "micro" ~ "Microcoleus",
                                    species == "nfix" ~ "Other N Fixers")) %>% 
  dplyr::mutate(model = case_when(model == "nitratefit" ~ "Nitrate",
                                    model == "phosfit" ~ "Phosphate",
                                    model == "ammoniumfit" ~ "Ammonium",
                                    model == "DINfit" ~ "DIN"))

clean.model.indices$model <- factor(  #Manually order model name 
  clean.model.indices$model,
  levels = c("Nitrate", "Phosphate", "Ammonium", "DIN")
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
ggplot(subset(clean.model.indices, metric %in% c("r2", "RMSE")), aes(x = value, y = model, shape = species, color = species)) +
  facet_wrap(~ metric, scales = "free_x") +
  geom_point(position = position_dodge(width = 0.6), #position_dodge seps species apart
             size = 3) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2,
                 position = position_dodge(width = 0.6)) +
  scale_y_discrete(limits = rev(levels(clean.model.indices$model)[1:4])) + #Reverses order of yaxis
  theme_bw() +
  colScale + shapScale +
  labs(x = "Metric Value",
       y = "Model Name",
       title = "River-Wide Nutrient Isolation Goodness-of-Fit",
       shape = "Species",
       color = "Species")


ggplot(subset(subset(clean.model.indices, category %in% "River-Wide"), metric == "RMSE"), aes(x = value, y = model, shape = species, color = species)) +
  geom_point(position = position_dodge(width = 0.6), #position_dodge seps species apart
             size = 3) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2,
                 position = position_dodge(width = 0.6)) +
  scale_y_discrete(limits = rev(levels(clean.model.indices$model)[1:4])) +
  colScale + shapScale +
  coord_cartesian(xlim = c(0, 5)) +
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
