#Riverwide Model Fit Calculations
library(ggplot2)
library(ggpubr)
library(tidyverse)

#OBSERVED DATA VS POSTERIORS
###############------------------------------------------------------------------
#Extract the relevant model here. fit.m4 = all variables,
#fit.m5 = biotic interactions, and fit.m6 = abiotic effects.
posteriors <- rstan::extract(fit.m4)[["n"]] #array indexed by iterations, species #, time

###############------------------------------------------------------------------
#Extract observed data vectors from alltaxatime, raw data object
y_micro <- alltaxatime$microcoleus
y_ana <- alltaxatime$anabaena_cylindrospermum
y_green <- alltaxatime$green_algae
y_nfix <- alltaxatime$other_nfixers

#Put them into a single matrix
y_obs_river <- rbind(y_green, y_micro, y_ana, y_nfix) #Dimensions: Species, time

###############------------------------------------------------------------------
#Begin calculating fit indices per species
iter <- dim(posteriors)[1]  # Number of iterations
species <- dim(posteriors)[2]  # Number of species
time <- dim(posteriors)[3]  # Time steps

R2 <- matrix(NA, iter, species) #Create empty matrix for R2 values per iteration, per species

for (s in 1:species) {
  
  y_obs <- y_obs_river[s, ]
  obs_index <- which(y_obs != -99) #Remove weeks where we did not collect field data
  
  for (i in 1:iter) {
    
    yhat <- posteriors[i, s, obs_index] #In the time index, take out weeks with no field data
    yobs <- y_obs[obs_index]
    
    var_fit <- var(yhat)
    var_res <- var(yobs - yhat)
    
    R2[i, s] <- var_fit / (var_fit + var_res)
  }
}

# Posterior summarize 
mean_R2 <- apply(R2, 2, mean) #2 stands for applying function over the columms
CI_R2 <- apply(R2, 2, quantile, c(0.025, 0.975))

mean_R2
CI_R2


#POSTERIORS VS PREDICTED
###############------------------------------------------------------------------
#Read in simulated data. But you must check that the simulation ran the same model
#(fit.m4, m5, m6) as the one being analyzed here
source(here::here("data_analysis/Predictions.R"))

#Comparing fit.m4: All variables included
S <- dim(latent)[1]
G <- dim(latent)[2]
T <- dim(latent)[3]

RMSE <- matrix(NA, S, G)
R2 <- matrix(NA, S, G)

for (g in 1:G) {
  for (s in 1:S) {
    
    z_lat <- latent[s, g, ]
    z_pred <- pred[s, g, ]
    
    SS_res <- sum((z_lat - z_pred)^2)
    SS_tot <- sum((z_lat - mean(z_lat))^2)
    
    RMSE[s, g] <- sqrt(mean((z_lat - z_pred)^2))
    R2[s, g] <- 1 - SS_res / SS_tot
  }
}

#Summarize RMSE
apply(RMSE, 2, mean)
apply(RMSE, 2, quantile, c(0.025, 0.975))

#Summarise R2
apply(R2, 2, mean)
apply(R2, 2, quantile, c(0.025, 0.975))





