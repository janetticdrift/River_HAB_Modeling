###########################
#Making Predictions: Running Simulations Using Latent States
###########################
#


#Packages for tidying and visualizing prediction simulations
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(abind)

#Read in cleaned latent states, params2_all. This file also reads in the cleaning_HAB.R file
source(here::here("data_analysis/Compare Obs Vs Modeled Outputs/River_Wide_isolation_model_vs_real .R"))

#Read in latent states and effect coefficients
nitratelatent <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_Nitrate_predictions.rds"))
phosphatelatent <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_Phos_predictions.rds"))
ammoniumlatent <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_Ammonium_predictions.rds"))
DINlatent <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_DIN_predictions.rds"))

#----------------------------------------------------------------------------
#River-wide: Nitrate only

#Pull out community abundances and demographics 
x <- nitratelatent
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- x[["Beta"]][,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(runs, 4, time)) #iterations, species, time

#Pull out environmental effects
Ntheta <- x[["Ntheta"]][,]
nitrate <- stand_nut$nitrate_mg_N_L[1:time]

Dtheta <- x[["Dtheta"]][,]
dis <- discharge$stand_discharge[1:time]

Ttheta <- x[["Ttheta"]][,]
temp <- stand_nut$temp_C[1:time]

Ctheta <- x[["Ctheta"]][,]
cond <- stand_nut$cond_uS_cm[1:time]

Rtheta <- x[["Rtheta"]][,]
rad <- swradiation$stand_rad[1:time]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  nTheta <- Ntheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          nTheta[s]*nitrate[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
    }
  }
}

#Create datafram for model checking
modelcheck_2022 <- n

#Create dataframe
sims2022median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022median, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))


#
#Historical Predictions, 2023
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,14] #iterations, time, species #

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time))

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[14:(13+time)]

dis <- discharge$stand_discharge[14:(13+time)]

temp <- stand_nut$temp_C[14:(13+time)]

cond <- stand_nut$cond_uS_cm[14:(13+time)]

rad <- swradiation$stand_rad[14:(13+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  nTheta <- Ntheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          nTheta[s]*nitrate[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
      
    }
  }
}

#Create datafram for model checking
modelcheck_2023 <- n

#Create dataframe
sims2023median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2023lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2023 <- left_join(sims2023median, sims2023lquant, by=c("Species", "time")) %>%
  left_join(., sims2023uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 24, year = 2023) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))


#
#Historical Predictions, 2024
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,29] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 17 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time))

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[29:(28+time)]

dis <- discharge$stand_discharge[29:(28+time)]

temp <- stand_nut$temp_C[29:(28+time)]

cond <- stand_nut$cond_uS_cm[29:(28+time)]

rad <- swradiation$stand_rad[29:(28+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  nTheta <- Ntheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients too
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          nTheta[s]*nitrate[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
    }
  }
}

#Create datafram for model checking
modelcheck_2024 <- n

#Create dataframe
sims2024median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2024 <- left_join(sims2024median, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 24, year = 2024) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



#Join together simulation data
simsnitrate <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(median = Abundance)

#Compile model check dataframes into a single full timeseries matrix
predictives.river.nitrate <- abind(exp(modelcheck_2022), exp(modelcheck_2023), 
                        exp(modelcheck_2024), along = 3)

#Save predictive output of All Variables model
saveRDS(predictives.river.nitrate, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Nitrate.rds"))



#----------------------------------------------------------------------------
#River-wide: Phosphate Variable

#Pull out community abundances and demographics 
x <- phosphatelatent
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- x[["Beta"]][,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(runs, 4, time)) #iterations, species, time

#Pull out environmental effects
Ptheta <- x[["Ptheta"]][,]
phos <- stand_nut$oPhos_ug_P_L[1:time]

Dtheta <- x[["Dtheta"]][,]
dis <- discharge$stand_discharge[1:time]

Ttheta <- x[["Ttheta"]][,]
temp <- stand_nut$temp_C[1:time]

Ctheta <- x[["Ctheta"]][,]
cond <- stand_nut$cond_uS_cm[1:time]

Rtheta <- x[["Rtheta"]][,]
rad <- swradiation$stand_rad[1:time]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  pTheta <- Ptheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          pTheta[s]*phos[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
    }
  }
}

#Create datafram for model checking
modelcheck_2022 <- n

#Create dataframe
sims2022median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022median, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))


#
#Historical Predictions, 2023
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,14] #iterations, time, species #

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time))

#Pull out environmental effects
phos <- stand_nut$oPhos_ug_P_L[14:(13+time)]

dis <- discharge$stand_discharge[14:(13+time)]

temp <- stand_nut$temp_C[14:(13+time)]

cond <- stand_nut$cond_uS_cm[14:(13+time)]

rad <- swradiation$stand_rad[14:(13+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  pTheta <- Ptheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          pTheta[s]*phos[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
      
    }
  }
}

#Create datafram for model checking
modelcheck_2023 <- n

#Create dataframe
sims2023median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2023lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2023 <- left_join(sims2023median, sims2023lquant, by=c("Species", "time")) %>%
  left_join(., sims2023uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 24, year = 2023) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))


#
#Historical Predictions, 2024
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,29] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 17 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time))

#Pull out environmental effects
phos <- stand_nut$oPhos_ug_P_L[29:(28+time)]

dis <- discharge$stand_discharge[29:(28+time)]

temp <- stand_nut$temp_C[29:(28+time)]

cond <- stand_nut$cond_uS_cm[29:(28+time)]

rad <- swradiation$stand_rad[29:(28+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  pTheta <- Ptheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients too
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          pTheta[s]*phos[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
    }
  }
}

#Create datafram for model checking
modelcheck_2024 <- n

#Create dataframe
sims2024median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2024 <- left_join(sims2024median, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 24, year = 2024) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



#Join together simulation data
simsphosphate <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(median = Abundance)

#Compile model check dataframes into a single full timeseries matrix
predictives.river.phosphate <- abind(exp(modelcheck_2022), exp(modelcheck_2023), 
                                   exp(modelcheck_2024), along = 3)

#Save predictive output of All Variables model
saveRDS(predictives.river.phosphate, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Phosphate.rds"))


#-----------------------------------------------------------------------------
#River-wide: Ammonium Variable

#Pull out community abundances and demographics 
x <- ammoniumlatent
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- x[["Beta"]][,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(runs, 4, time)) #iterations, species, time

#Pull out environmental effects
Atheta <- x[["Atheta"]][,]
ammonium <- stand_nut$ammonium_mg_N_L[1:time]

Dtheta <- x[["Dtheta"]][,]
dis <- discharge$stand_discharge[1:time]

Ttheta <- x[["Ttheta"]][,]
temp <- stand_nut$temp_C[1:time]

Ctheta <- x[["Ctheta"]][,]
cond <- stand_nut$cond_uS_cm[1:time]

Rtheta <- x[["Rtheta"]][,]
rad <- swradiation$stand_rad[1:time]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          aTheta[s]*ammonium[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
    }
  }
}

#Create datafram for model checking
modelcheck_2022 <- n

#Create dataframe
sims2022median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022median, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))


#
#Historical Predictions, 2023
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,14] #iterations, time, species #

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time))

#Pull out environmental effects
ammonium <- stand_nut$ammonium_mg_N_L[14:(13+time)]

dis <- discharge$stand_discharge[14:(13+time)]

temp <- stand_nut$temp_C[14:(13+time)]

cond <- stand_nut$cond_uS_cm[14:(13+time)]

rad <- swradiation$stand_rad[14:(13+time)]

for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          aTheta[s]*ammonium[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
      
    }
  }
}

#Create datafram for model checking
modelcheck_2023 <- n

#Create dataframe
sims2023median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2023lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2023 <- left_join(sims2023median, sims2023lquant, by=c("Species", "time")) %>%
  left_join(., sims2023uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 24, year = 2023) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))


#
#Historical Predictions, 2024
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,29] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 17 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time))

#Pull out environmental effects
ammonium <- stand_nut$ammonium_mg_N_L[29:(28+time)]

dis <- discharge$stand_discharge[29:(28+time)]

temp <- stand_nut$temp_C[29:(28+time)]

cond <- stand_nut$cond_uS_cm[29:(28+time)]

rad <- swradiation$stand_rad[29:(28+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients too
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          aTheta[s]*ammonium[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
    }
  }
}

#Create datafram for model checking
modelcheck_2024 <- n

#Create dataframe
sims2024median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2024 <- left_join(sims2024median, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 24, year = 2024) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



#Join together simulation data
simsammonium <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(median = Abundance)

# Plot 2: Only show Anabaena + Microcoleus
ggplot(simsnitrate, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsnitrate,
                               median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsnitrate,
                             median = ifelse(Species %in% c("Anabaena", "Microcoleus"), 
                                             median, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 5)) +
  coord_cartesian(ylim = c(0,16)) +
  labs(x = "Date", y = "Percent Cover (%)") +
  colScale + filScale + linScale + theme_bw()

#Compile model check dataframes into a single full timeseries matrix
predictives.river.ammonium <- abind(exp(modelcheck_2022), exp(modelcheck_2023), 
                                     exp(modelcheck_2024), along = 3)

#Save predictive output of All Variables model
saveRDS(predictives.river.ammonium, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Ammonium.rds"))

#-----------------------------------------------------------------------------
#River-wide: DIN Variable

#Pull out community abundances and demographics 
x <- DINlatent
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- x[["Beta"]][,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(runs, 4, time)) #iterations, species, time

#Pull out environmental effects
DINtheta <- x[["DINtheta"]][,]
DIN <- stand_nut$DIN[1:time]

Dtheta <- x[["Dtheta"]][,]
dis <- discharge$stand_discharge[1:time]

Ttheta <- x[["Ttheta"]][,]
temp <- stand_nut$temp_C[1:time]

Ctheta <- x[["Ctheta"]][,]
cond <- stand_nut$cond_uS_cm[1:time]

Rtheta <- x[["Rtheta"]][,]
rad <- swradiation$stand_rad[1:time]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  dinTheta <- DINtheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          dinTheta[s]*DIN[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
    }
  }
}

#Create datafram for model checking
modelcheck_2022 <- n

#Create dataframe
sims2022median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022median, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))


#
#Historical Predictions, 2023
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,14] #iterations, time, species #

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time))

#Pull out environmental effects
DIN <- stand_nut$DIN[14:(13+time)]

dis <- discharge$stand_discharge[14:(13+time)]

temp <- stand_nut$temp_C[14:(13+time)]

cond <- stand_nut$cond_uS_cm[14:(13+time)]

rad <- swradiation$stand_rad[14:(13+time)]

for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  dinTheta <- DINtheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          dinTheta[s]*DIN[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
      
    }
  }
}

#Create datafram for model checking
modelcheck_2023 <- n

#Create dataframe
sims2023median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2023lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2023 <- left_join(sims2023median, sims2023lquant, by=c("Species", "time")) %>%
  left_join(., sims2023uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 24, year = 2023) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))


#
#Historical Predictions, 2024
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,29] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 17 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time))

#Pull out environmental effects
DIN <- stand_nut$DIN[29:(28+time)]

dis <- discharge$stand_discharge[29:(28+time)]

temp <- stand_nut$temp_C[29:(28+time)]

cond <- stand_nut$cond_uS_cm[29:(28+time)]

rad <- swradiation$stand_rad[29:(28+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  #Pull env covariates
  dinTheta <- DINtheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients too
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
                          dinTheta[s]*DIN[t-1] +
                          dTheta[s]*dis[t-1] + tTheta[s]*temp[t-1] + 
                          cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1], 
                        sd = sigma[s])
    }
  }
}

#Create datafram for model checking
modelcheck_2024 <- n

#Create dataframe
sims2024median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2024), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2024 <- left_join(sims2024median, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 24, year = 2024) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



#Join together simulation data
simsDIN <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(median = Abundance)

#Compile model check dataframes into a single full timeseries matrix
predictives.river.DIN <- abind(exp(modelcheck_2022), exp(modelcheck_2023), 
                                    exp(modelcheck_2024), along = 3)

#Save predictive output of All Variables model
saveRDS(predictives.river.DIN, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_DIN.rds"))
