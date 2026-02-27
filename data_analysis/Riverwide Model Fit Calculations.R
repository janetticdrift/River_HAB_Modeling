#Riverwide Model Fit Calculations, 2022
library(ggplot2)
library(ggpubr)
library(tidyverse)

#RAW DATA VS POSTERIORS
###############------------------------------------------------------------------
#Extract the relevant model here. fit.m4 = all variables,
#fit.m5 = biotic interactions, and fit.m6 = abiotic effects.
posteriors <- rstan::extract(fit.m4)[["n"]] #array indexed by iterations, species #, time

#Pull out each species result. 
#1 = green algae, 2 = microcoleus, 3 = anabaena, 4 = other nitrogen fixers
y_rep_micro <- posteriors[,2,] 
y_rep_ana <- posteriors[,3,] 
y_rep_green <- posteriors[,1,] 
y_rep_nfix <- posteriors[,4,] 

###############------------------------------------------------------------------
#Extract observed data vectors from alltaxatime, raw data object
y_micro <- alltaxatime$microcoleus
y_ana <- alltaxatime$anabaena_cylindrospermum
y_green <- alltaxatime$green_algae
y_nfix <- alltaxatime$other_nfixers

###############------------------------------------------------------------------
#Begin calculating fit index per model iteration.
#Plug in current species of interest
y_obs <- y_nfix # y_observed data
y_rep <- y_rep_nfix # y_replicated data

obs_index <- which(y_obs != -99) # Remove placeholder values from comparison
Iter <- nrow(y_rep) # number of iterations in the model
R2 <- numeric(Iter) # Create empty vector to be filled

#Calculate R^2 for each species
for (i in 1:Iter) {
  yhat <- y_rep[i, obs_index] # Latent states
  yobs <- y_obs[obs_index] # Observed states
  
  var_fit <- var(yhat) #Calculate variance of latent states
  var_res <- var(yobs - yhat) #calculate residual variance (diff between obs and predicted)
  
  R2[i] <- var_fit / (var_fit + var_res)
}

# Posterior summary
mean_R2 <- mean(R2)
CI_R2 <- quantile(R2, c(0.025, 0.975))

mean_R2
CI_R2
