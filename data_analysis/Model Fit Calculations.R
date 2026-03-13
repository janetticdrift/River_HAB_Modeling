#Model Fit Calculations
library(ggplot2)
library(ggpubr)
library(tidyverse)

#Must read in alltaxatime
#Must read in atx 

#RIVERWIDE

#OBSERVED DATA VS POSTERIORS
###############------------------------------------------------------------------
#Extract observed data vectors from alltaxatime, raw data object
y_micro <- alltaxatime$microcoleus
y_ana <- alltaxatime$anabaena_cylindrospermum
y_green <- alltaxatime$green_algae
y_nfix <- alltaxatime$other_nfixers

#Put them into a single matrix
y_obs_river <- rbind(y_green, y_micro, y_ana, y_nfix) #Dimensions: Species, time

###############------------------------------------------------------------------
#Read in modeled data
allfit <- readRDS(here::here("data/Riverwide_AllVar_predictions.rds"))
bioticfit <- readRDS(here::here("data/Riverwide_Biotic_predictions.rds"))
abioticfit <- readRDS(here::here("data/Riverwide_Abiotic_predictions.rds"))
abioticnonutfit <- readRDS(here::here("data/Riverwide_AbioticNonut_predictions.rds"))

#Extract the relevant model here. fit.m4 = all variables,
#fit.m5 = biotic interactions, and fit.m6 = abiotic effects.
posteriors <- TMfit[["n"]] #array indexed by iterations, species #, time


###############------------------------------------------------------------------
#BAYESIAN R2 CALCULATIONS

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


#####################################################################################
#####################################################################################
#####################################################################################

#WITHIN MAT

#Read in data
TMfit <- readRDS(here::here("data/WithinMat_Micro_predictions.rds"))
TAfit <- readRDS(here::here("data/WithinMat_Ana_predictions.rds"))

#OBSERVED DATA VS POSTERIORS
###############------------------------------------------------------------------
#Extract the relevant model here. fit.m1 = microcoleus, fit.m2 = anabaena.
posteriors_mat <- TMfit[["n"]] #array indexed by iterations, species #, time

###############------------------------------------------------------------------
#Extract observed data vectors from alltaxatime, raw data object
y_ana <- matalltaxaM$Anabaena                  #1
y_epi <- matalltaxaM$`Epithemia Diatoms`       #2
y_geit <- matalltaxaM$Geitlerinema             #3
y_green <- matalltaxaM$`Green Algae`           #4
y_micro <- matalltaxaM$Microcoleus             #5
y_nonepi <- matalltaxaM$`Non-Epithemia Diatoms`#6
y_nos <- matalltaxaM$Nostoc                    #7
y_cocc <- matalltaxaM$`Other Coccoids`         #8
y_rare <- matalltaxaM$Rare                     #9

#Put them into a single matrix
y_obs_mat <- rbind(y_ana, y_epi, y_geit, y_green, y_micro,
                   y_nonepi, y_nos, y_cocc, y_rare) #Dimensions: Species, time


###############------------------------------------------------------------------
#Bayesian R2 comparing latent vs observed states
iter <- dim(posteriors_mat)[1]  # Number of iterations
species <- dim(posteriors_mat)[2]  # Number of species
time <- dim(posteriors_mat)[3]  # Time steps

names<- c("Anabaena", "Epithemia", "Geitlerinema", 
          "Green Algae", "Microcoleus",
          "Non-Epithemia", "Nostoc", "Other Coccoids",
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
apply(R2, 2, median)[c(1:3)] #2 stands for applying function over the columms
#Credible Interval
apply(R2, 2, quantile, c(0.025, 0.975))


#POSTERIORS VS PREDICTED
###############------------------------------------------------------------------
#Read in simulated data. Already back-transformed 
TM.pred <- readRDS(here::here("data/WithinMat_Pred_TM.rds"))

#Set current predictive data
predictives.mats <- TM.pred

#Back-transform latent abundance data
logposteriors.mats <- exp(posteriors_mat)


#Comparing fit.m1: All variables included
iter <- dim(logposteriors.mats)[1]  # Number of iterations 
species <- dim(logposteriors.mats)[2]  # Number of species 
time <- dim(logposteriors.mats)[3]  # Time steps

names<- c("Anabaena", "Epithemia", "Geitlerinema", 
          "Green Algae", "Microcoleus",
          "Non-Epithemia", "Nostoc", "Other Coccoids",
          "Rare")

#Create empty matrices for storing fit index values
RMSE <- matrix(NA, iter, species, dimnames = list(NULL, names))
R2 <- matrix(NA, iter, species, dimnames = list(NULL, names))

for (s in 1:species) {
  for (i in 1:iter) {
    
    y <- logposteriors.mats[i, s, ]
    y_pred <- predictives.mats[i, s, ]
    
    RMSE[i, s] <- sqrt(mean((y - y_pred)^2)) #Calculate RMSE per species iteration
    
    R2[i, s] <- cor(y, y_pred)^2 #Calculate R2 per species iteration
  }
}

#Summarize RMSE
apply(RMSE, 2, mean)
apply(RMSE, 2, quantile, c(0.025, 0.975))

#Summarise R2
apply(R2, 2, mean)
apply(R2, 2, quantile, c(0.025, 0.975))



#Correlation between latent Anabaena and Epithemia
cors <- apply(logposteriors.mats[, c(1,2), ], 1, function(x){
  cor(x[1,], x[2,])
})
mean(cors)

#Correlation between predicted Anabaena and Epithemia
cors <- apply(predictives.mats[, c(1,2), ], 1, function(x){
  cor(x[1,], x[2,])
})
mean(cors)

# sync_fun <- function(mat) {
#   # mat = species × time
#   total <- colSums(mat)
#   var(total) / (sum(apply(mat, 1, sd))^2)
# }
# 
# sync <- apply(logposteriors.mats[, c(1,2), 1:13], 1, function(x) sync_fun(x))
# mean(sync)

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

