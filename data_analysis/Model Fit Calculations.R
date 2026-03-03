#Model Fit Calculations
library(ggplot2)
library(ggpubr)
library(tidyverse)

#RIVERWIDE

#OBSERVED DATA VS POSTERIORS
###############------------------------------------------------------------------
#Extract the relevant model here. fit.m4 = all variables,
#fit.m5 = biotic interactions, and fit.m6 = abiotic effects.
posteriors <- rstan::extract(fit.m6)[["n"]] #array indexed by iterations, species #, time

###############------------------------------------------------------------------
#Extract observed data vectors from alltaxatime, raw data object
y_micro <- alltaxatime$microcoleus
y_ana <- alltaxatime$anabaena_cylindrospermum
y_green <- alltaxatime$green_algae
y_nfix <- alltaxatime$other_nfixers

#Put them into a single matrix
y_obs_river <- rbind(y_green, y_micro, y_ana, y_nfix) #Dimensions: Species, time

###############------------------------------------------------------------------
#Bayesian R2 comparing latent vs observed states
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
  #Mean R2
apply(R2, 2, mean) #2 stands for applying function over the columms
  #Credible Interval
apply(R2, 2, quantile, c(0.025, 0.975))


#POSTERIORS VS PREDICTED
###############------------------------------------------------------------------
#Read in simulated data. But you must check that the simulation ran the same model
#(fit.m4, m5, m6) as the one being analyzed here
source(here::here("data_analysis/Predictions.R"))
#Relevant dataframe is "predictives"

#Comparing fit.m4: All variables included
iter <- dim(posteriors)[1]  # Number of iterations 
species <- dim(posteriors)[2]  # Number of species 
time <- dim(posteriors)[3]  # Time steps

#Create empty matrices for storing fit index values
RMSE <- matrix(NA, iter, species)
R2 <- matrix(NA, iter, species)

for (s in 1:species) {
  for (i in 1:iter) {
    
    y <- posteriors[i, s, ]
    y_pred <- predictives[i, s, ]
    
    RMSE[i, s] <- sqrt(mean((y - y_pred)^2)) #Calculate RMSE per species iteration
    
    SS_res <- sum((y - y_pred)^2)
    SS_tot <- sum((y - mean(y))^2)
    
    R2[i, s] <- 1 - SS_res / SS_tot #Calculate R2 per species iteration
  }
}

#Summarize RMSE
apply(RMSE, 2, mean)
apply(RMSE, 2, quantile, c(0.025, 0.975))

#Summarise R2
apply(R2, 2, mean)
apply(R2, 2, quantile, c(0.025, 0.975))


#####################################################################################
#####################################################################################
#####################################################################################

#WITHIN MAT
#OBSERVED DATA VS POSTERIORS
###############------------------------------------------------------------------
#Extract the relevant model here. fit.m1 = microcoleus, fit.m2 = anabaena.
posteriors_mat <- rstan::extract(fit.m1)[["n"]] #array indexed by iterations, species #, time

###############------------------------------------------------------------------
#Extract observed data vectors from alltaxatime, raw data object
y_ana <- matalltaxaM$Anabaena                  #1
y_epi <- matalltaxaM$`Epithemia Diatoms`       #2
y_geit <- matalltaxaM$Geitlerinema             #3
y_green <- matalltaxaM$`Green Algae`           #4
y_lep <- matalltaxaM$Leptolyngbya              #5
y_micro <- matalltaxaM$Microcoleus             #6
y_nonepi <- matalltaxaM$`Non-Epithemia Diatoms`#7
y_nos <- matalltaxaM$Nostoc                    #8
y_osc <- matalltaxaM$Oscillatoria              #9
y_cocc <- matalltaxaM$`Other Coccoids`         #10
y_rare <- matalltaxaM$Rare                     #11

#Put them into a single matrix
y_obs_mat <- rbind(y_ana, y_epi, y_geit, y_green, y_lep, y_micro,
                   y_nonepi, y_nos, y_osc, y_cocc, y_rare) #Dimensions: Species, time


###############------------------------------------------------------------------
#Bayesian R2 comparing latent vs observed states
iter <- dim(posteriors_mat)[1]  # Number of iterations
species <- dim(posteriors_mat)[2]  # Number of species
time <- dim(posteriors_mat)[3]  # Time steps

names<- c("Anabaena", "Epithemia", "Geitlerinema", 
          "Green Algae", "Leptolyngbya", "Microcoleus",
          "Non-Epithemia", "Nostoc", "Oscillatoria", "Other Coccoids",
          "Rare")

R2 <- matrix(NA, iter, species, dimnames = list(NULL, names)) #Create empty matrix for R2 values per iteration, per species

for (s in 1:species) {
  
  y_obs <- y_obs_mat[s, ]
  obs_index <- which(y_obs != -99) #Remove weeks where we did not collect field data
  
  for (i in 1:iter) {
    
    yhat <- posteriors_mat[i, s, obs_index] #In the time index, take out weeks with no field data
    yobs <- y_obs[obs_index]
    
    var_fit <- var(yhat)
    var_res <- var(yobs - yhat)
    
    R2[i, s] <- var_fit / (var_fit + var_res)
  }
}

# Posterior summarize 
#Mean R2
apply(R2, 2, mean)[c(4,6,7)] #2 stands for applying function over the columms
#Credible Interval
apply(R2, 2, quantile, c(0.025, 0.975))


#POSTERIORS VS PREDICTED
###############------------------------------------------------------------------
#Read in simulated data. But you must check that the simulation ran the same model
#(fit.m4, m5, m6) as the one being analyzed here
source(here::here("data_analysis/Predictions.R"))
#Relevant dataframe is "predictives"

#Comparing fit.m1: All variables included
iter <- dim(posteriors_mat)[1]  # Number of iterations 
species <- dim(posteriors_mat)[2]  # Number of species 
time <- dim(posteriors_mat)[3]  # Time steps

names<- c("Anabaena", "Epithemia", "Geitlerinema", 
          "Green Algae", "Leptolyngbya", "Microcoleus",
          "Non-Epithemia", "Nostoc", "Oscillatoria", "Other Coccoids",
          "Rare")

#Create empty matrices for storing fit index values
RMSE <- matrix(NA, iter, species, dimnames = list(NULL, names))
R2 <- matrix(NA, iter, species, dimnames = list(NULL, names))

for (s in 1:species) {
  for (i in 1:iter) {
    
    y <- posteriors_mat[i, s, ]
    y_pred <- predictives_mats[i, s, ]
    
    RMSE[i, s] <- sqrt(mean((y - y_pred)^2)) #Calculate RMSE per species iteration
    
    SS_res <- sum((y - y_pred)^2)
    SS_tot <- sum((y - mean(y))^2)
    
    R2[i, s] <- 1 - SS_res / SS_tot #Calculate R2 per species iteration
  }
}

#Summarize RMSE
apply(RMSE, 2, mean)
apply(RMSE, 2, quantile, c(0.025, 0.975))

#Summarise R2
apply(R2, 2, mean)
apply(R2, 2, quantile, c(0.025, 0.975))
