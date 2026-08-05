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

                                #Toxins
#Read in OBSERVED states of Target Microcoleus toxin concentrations
y_obs_tox_TM <- readRDS(here::here("data/Outputs for Sims and Model Fits/obs_toxins_data_TM.rds"))
y_obs_tox_TAC <- readRDS(here::here("data/Outputs for Sims and Model Fits/obs_toxins_data_TAC.rds"))

#Read in LATENT states of Toxin Models
toxriverfitTM <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_predictions.rds"))
toxriverfitTAC <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TAC_River_predictions.rds"))
toxmatfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_Mat_predictions.rds"))

#Read in SIMULATED states of Toxin Models
pred.toxriverfitTM <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Toxins_TM_Pred_RiverWide.rds"))
pred.toxriverfitTAC <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Toxins_TAC_Pred_RiverWide.rds"))
pred.toxmatfit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Toxins_Pred_Withinmat.rds"))


###############------------------------------------------------------------------
#Create lists of the models to iterate through in for loop
river_species <- c("green", "micro", "ana", "nfix")
mat_species <- c("ana", "epi", "geit")
tox_congener <- "Anatoxin"

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
  ),
  list(
    model = "toxriverfitTM",
    category = "River-Wide",
    y_obs = y_obs_tox_TM$ATX_all_ug_afdm_g,
    post = toxriverfitTM[["tox_raw"]],
    pred = pred.toxriverfitTM,
    species = tox_congener
  ),
  list(
    model = "toxriverfitTAC",
    category = "River-Wide",
    y_obs = y_obs_tox_TAC$ATX_all_ug_afdm_g,
    post = toxriverfitTAC[["tox_raw"]],
    pred = pred.toxriverfitTAC,
    species = tox_congener
  ),
  list(
    model = "toxmatfit",
    category = "Within-Mat",
    y_obs = y_obs_tox_TM$ATX_all_ug_afdm_g,
    post = toxmatfit[["tox_raw"]],
    pred = pred.toxmatfit,
    species = tox_congener
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
  posteriors <- model$post          #Posterior values
  predictives <- model$pred         #Predicted values
  species_names <- model$species    #Names of species in this model
  
  #For toxins with a 2D posterior (iterations x time), convert to iterations x 1 x time
  if(length(dim(posteriors)) == 2){
    posteriors <- array(posteriors,
      dim = c(dim(posteriors)[1], 1, dim(posteriors)[2]) #Add in a third dimension in the "species" middle spot
    )
    predictives <- array(predictives,
      dim = c(dim(predictives)[1], 1, dim(predictives)[2])
    )
    
    y_obs <- matrix(y_obs, nrow = 1)  #Convert the y_obs vector into a matrix
  }
  
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
    obs_index <- which(y_obs_s != -99) #Do not include weeks where we did not observe field data
    
    for (i in 1:iter) {
      #Bayesian R2      
      yhat <- posteriors[i, s, obs_index] #In the time index position, take out weeks with no field observations
      yobs <- y_obs_s[obs_index]          #Take out weeks where we did not have field observations
      
      var_fit <- var(yhat)
      var_res <- var(yobs - yhat)
      
      BayesR2_vals[i, s] <- var_fit / (var_fit + var_res)
      
      #RMSE
      y <- posteriors[i, s, ]
      y_pred <- predictives[i, s, ]

      RMSE_vals[i, s] <- sqrt(mean((y - y_pred)^2)) #Calculate RMSE per species iteration

      range_y <- max(y) - min(y)                    #Calculate the range of values

      NRMSE_vals[i, s] <- RMSE_vals[i, s] / range_y #Use range to normalize RMSE

      #R2
      R2_vals[i, s] <- cor(y, y_pred)^2             #Calculate R2 per species iteration
      
      
    }
  } 
  
  #Summarize metrics across all iterations
  metrics <- list(
    "Bayesian R2" = BayesR2_vals,
    "RMSE" = RMSE_vals,
    "NRMSE" = NRMSE_vals,
    "r2" = R2_vals
  )
  
  #Calculate median and CIs for each metric
  for(j in names(metrics)){
    
    vals <- metrics[[j]]                            #Subset out values for the named metric
    
    med <- apply(vals, 2, mean)                   #Calculate mean of that metric
    ci <- apply(vals, 2, quantile, c(0.025, 0.975)) #Calculate CI of that metric
    
    for(s in 1:species){                            #Place mean and CI into the 'model.indices' dataframe
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
      
      counter <- counter + 1                      #Add to counter to run through the next model
      
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
                                    species == "nfix" ~ "Other N Fixers",
                                    species == "Anatoxin" ~ "Anatoxin")) %>% 
  dplyr::mutate(model = case_when(model == "allfit" ~ "Percent Cover: All Variables",
                                    model == "bioticfit" ~ "Percent Cover: Biotic Only",
                                    model == "abioticfit" ~ "Percent Cover: Abiotic Only",
                                    model == "abioticnonutfit" ~ "Percent Cover: Abiotic Minus Nutrients",
                                    model == "TMfit" ~ "Microscopy: All Variables",
                                  model == "toxriverfitTM" ~ "Microcoleus Mat Toxins",
                                  model == "toxriverfitTAC" ~ "Anabaena Mat Toxins",
                                  model == "toxmatfit" ~ "Microcoleus Mat Toxins")) %>% 
   dplyr::filter(species %in% c("Anabaena", "Microcoleus", "Anatoxin"))

clean.model.indices$model <- factor(  #Manually order model name 
  clean.model.indices$model,
  levels = c("Percent Cover: All Variables", "Percent Cover: Biotic Only", 
             "Percent Cover: Abiotic Only", 
             "Percent Cover: Abiotic Minus Nutrients", "Microscopy: All Variables",
             "Microcoleus Mat Toxins",
             "Anabaena Mat Toxins")
)
clean.model.indices$species <- factor(  #Manually order species name 
  clean.model.indices$species,
  levels = c("Anatoxin", "Anabaena", "Microcoleus")
)

#Create a color palette
mycols <- c(
  "Percent Cover: All Variables" = "#095dbb",
  "Percent Cover: Biotic Only" = "#00899a",
  "Percent Cover: Abiotic Only" = "#11AB7C",
  "Percent Cover: Abiotic Minus Nutrients" = "#37Ca61",
  "Microscopy: All Variables" = "#93C451",
  "Microcoleus Mat Toxins" = "#BE5103",
  "Anabaena Mat Toxins" = "#ff9913"
  # "Microcoleus Mat Toxins" = "#ffd5b6")
)
colScale <- scale_color_manual(
  name = "Model",
  values = mycols)

myshap <- c("River-Wide" = 16, "Within-Mat" = 17)
shapScale <- scale_shape_manual(
  name = "Data Source",
  values = myshap)

#Plot River-Wide Algae Metrics
ggplot(subset(clean.model.indices, metric %in% c("r2", "NRMSE")), 
       aes(x = value, y = species, shape = category, color = model)) +
  facet_wrap(~ metric, scales = "free_x") +
  geom_point(position = position_dodge(width = -0.6), #position_dodge separates species apart
             size = 3) +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2,
                 position = position_dodge(width = -0.6)) +
  colScale + shapScale +
  theme_bw() +
  labs(x = "Metric Value",
       y = "",
       title = "Goodness-of-Fit",
       shape = "Model",
       color = "Model")


# ggplot(subset(subset(clean.model.indices, category %in% "River-Wide"), metric == "Bayesian R2"), aes(x = value, y = species, shape = model, color = model)) +
#   geom_point(position = position_dodge(width = 0.6), #position_dodge seps species apart
#              size = 3) +
#   geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2,
#                  position = position_dodge(width = 0.6)) +
#   scale_y_discrete(limits = rev(levels(clean.model.indices$species)[1:4])) +
#   colScale + shapScale +
#   theme_bw() +
#   labs(title = "Bayesian R2",
#        x = "Metric value",
#        y = "Model",
#        shape = "Species",
#        color = "Species")
