#Model Fit Calculations
library(ggplot2)
library(ggpubr)
library(tidyverse)

#Must read in atx 

###############------------------------------------------------------------------
#Create empty dataframe for storing all calculated indices
model.indices <- data.frame(
  model = character(),
  category = character(),
  species = character(),
  metric = character(),
  value = numeric(),
  lwr = numeric(),
  upr = numeric(),
  stringsAsFactors = FALSE #Do not make categorical variables factors
)

###############------------------------------------------------------------------
#Read in all data needed for calculations

                               #River-Wide

#Read in OBSERVED states of River-Wide
obs_river_data <- readRDS(here::here("data/Model Fits/obs_river_data.rds"))
#Extract log-transformed observed data vectors from obs_river_data, 
y_micro <- obs_river_data$microcoleus
y_ana <- obs_river_data$anabaena_cylindrospermum
y_green <- obs_river_data$green_algae
y_nfix <- obs_river_data$other_nfixers
#Put them into a single matrix
y_obs_river <- rbind(y_green, y_micro, y_ana, y_nfix) #Dimensions: Species, time

#Read in LATENT states of River-Wide
allfit <- readRDS(here::here("data/Model Fits/Riverwide_AllVar_predictions.rds"))
bioticfit <- readRDS(here::here("data/Model Fits/Riverwide_Biotic_predictions.rds"))
abioticfit <- readRDS(here::here("data/Model Fits/Riverwide_Abiotic_predictions.rds"))
abioticnonutfit <- readRDS(here::here("data/Model Fits/Riverwide_AbioticNonut_predictions.rds"))

#Read in SIMULATED states of River-Wide
pred.allfit <- readRDS(here::here("data/Riverwide_Pred_AllVar.rds"))
pred.bioticfit <- readRDS(here::here("data/Riverwide_Pred_Biotic.rds"))
pred.abioticfit <- readRDS(here::here("data/Riverwide_Pred_Abiotic.rds"))
pred.abioticnonutfit <- readRDS(here::here("data/Riverwide_Pred_AbioticNoNut.rds"))


                               #Within-Mat

#Read in OBSERVED states of Within-Mat: TM
obs_mat_data <- readRDS(here::here("data/Model Fits/obs_TM_data.rds"))
#Extract log-transformed observed data vectors from obs_mat_data, raw data object
y_ana <- obs_mat_data$Anabaena                  #1
y_epi <- obs_mat_data$`Epithemia Diatoms`       #2
y_geit <- obs_mat_data$Geitlerinema             #3
#Put them into a single matrix
y_obs_mat <- rbind(y_ana, y_epi, y_geit) #Dimensions: Species, time

#Read in LATENT states of River-Wide: TM
TMfit <- readRDS(here::here("data/Model Fits/WithinMat_Micro_predictions.rds"))

#Read in SIMULATED states of River-Wide
TM.pred <- readRDS(here::here("data/WithinMat_Pred_TM.rds"))

#Create lists of models to iterate through
river_species <- c("green", "micro", "ana", "nfix")
mat_species <- c("ana", "epi", "geit")

model_list <- list(
  list(
    model = "allfit",          #The model name
    category = "River-Wide",  #Whether it used percent cover or microscopy data
    y_obs = y_obs_river,      #Reading in observed values
    post = allfit[["n"]],     #Reading in latent state posteriors
    pred = pred.allfit,       #Reading in simulated/predicted values
    species = river_species    #List of species names used in this model
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

counter <- 1

for (m in 1:length(model_list)) {
  
  model <- model_list[[m]]
  
  y_obs <- model$y_obs
  posteriors <- model$post          #For comparing obs with latent on same scale
  logposteriors <- exp(posteriors)  #For comparing latent with pred on same scale
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
      
      y <- logposteriors[i, s, ]
      y_pred <- predictives[i, s, ]
      
      RMSE_vals[i, s] <- sqrt(mean((y - y_pred)^2)) #Calculate RMSE per species iteration
 
                                         #R2      
      R2_vals[i, s] <- cor(y, y_pred)^2 #Calculate R2 per species iteration
      
    }
  } 
                                   #Calculate metrics
    #Consolidate metrics per iter
    metrics <- list(
      "Bayesian R2" = BayesR2_vals,
      "RMSE" = RMSE_vals, 
      "R2" = R2_vals
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

#####################################################################################
#####################################################################################
#####################################################################################

#TOXINS

#OBSERVED DATA VS POSTERIORS
###############------------------------------------------------------------------
#Extract observed data vectors from atx, raw data object
y_obs_atx <- anatoxin_data$ATX_all_ug_g


###############------------------------------------------------------------------
#Read in modeled data
atxfit <- readRDS(here::here("data/Anatoxin_AllVar_predictions.rds"))

#Extract the relevant model here. fit.m4 = all variables,
#fit.m5 = biotic interactions, and fit.m6 = abiotic effects.
posteriors <- atxfit[["tox"]] #array indexed by iterations, species #, time


###############------------------------------------------------------------------
#BAYESIAN R2 CALCULATIONS

#Bayesian R2 comparing latent vs observed states
iter <- dim(posteriors)[1]  # Number of iterations

R2 <- rep(NA, times = iter)  #Create empty matrix for R2 values per iteration, per species

  for (i in 1:iter) {
    
    y_obs <- y_obs_atx[ ]
    obs_index <- which(y_obs != -99) #Remove weeks where we did not collect field data
    
    yhat <- posteriors[i, obs_index] #In the time index, take out weeks with no field data
    yobs <- y_obs[obs_index]
    
    var_fit <- var(yhat)
    var_res <- var(yobs - yhat)
    
    R2[i] <- var_fit / (var_fit + var_res)
  }


# Posterior summarize 
#Mean R2
mean(R2) #2 stands for applying function over the columms
#Credible Interval
quantile(R2, c(0.025, 0.975))


#POSTERIORS VS PREDICTED
###############------------------------------------------------------------------
#LATENT VS PREDICTION RMSE & R2 CALCULATIONS

#Read in simulated data
pred.all.fit <- readRDS(here::here("data/Riverwide_Pred_AllVar.rds"))
pred.biotic.fit <- readRDS(here::here("data/Riverwide_Pred_Biotic.rds"))
pred.abiotic.fit <- readRDS(here::here("data/Riverwide_Pred_Abiotic.rds"))
pred.abioticnonut.fit <- readRDS(here::here("data/Riverwide_Pred_AbioticNoNut.rds"))

#Set current predictive data
predictives <- pred.abioticnonut.fit

#Back-transform latent abundance data
logposteriors <- exp(posteriors)


#Comparing All variables included
iter <- dim(logposteriors)[1]  # Number of iterations 
species <- dim(logposteriors)[2]  # Number of species 
time <- dim(logposteriors)[3]  # Time steps

#Create empty matrices for storing fit index values
RMSE <- matrix(NA, iter, species)
R2 <- matrix(NA, iter, species)

for (s in 1:species) {
  for (i in 1:iter) {
    
    y <- logposteriors[i, s, ]
    y_pred <- predictives[i, s, ]
    
    RMSE[i, s] <- sqrt(mean((y - y_pred)^2)) #Calculate RMSE per species iteration
    
    R2[i, s] <- cor(y, y_pred)^2 #Calculate R2 per species iteration
  }
}

#Summarize RMSE
apply(RMSE, 2, median)
apply(RMSE, 2, quantile, c(0.025, 0.975))

#Summarise R2
apply(R2, 2, median)
apply(R2, 2, quantile, c(0.025, 0.975))

