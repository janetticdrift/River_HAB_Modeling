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
source(here::here("data_analysis/Compare Obs Vs Modeled Outputs/River_Wide_isolation_model_vs_real.R"))

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
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022median, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 25, year = 2022) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
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

#Plot simulated predictions against latent states
# Plot 1: Only show Green Algae + Other N Fixers
p1 <- ggplot(simsnitrate, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsnitrate,
                               median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
             data = transform(simsnitrate,
                              median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), 
                                            median, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,72)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Nitrate Isolation") +
  colScale + filScale + theme_bw()

# Plot 2: Only show Anabaena + Microcoleus
p2 <- ggplot(simsnitrate, aes(x = model_date, y = median)) +
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
  colScale + filScale + theme_bw()
  

# Combine plots and collect legends
  # See colScale code on line 22 to add/remove taxa from the Taxa list, using the breaks function
(p1 / p2) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

# Combine Microcoleus/Anabaena plot with environmental variables plot
(p2 / envplot) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical")


#Compile model check dataframes into a single full timeseries matrix
predictives.river.nitrate <- abind(exp(modelcheck_2022), exp(modelcheck_2023), 
                        exp(modelcheck_2024), along = 3)

#Save predictive output of All Variables model
saveRDS(predictives.river.nitrate, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Nitrate.rds"))



#----------------------------------------------------------------------------
#River-wide: Biotic Variables

#Pull out community abundances and demographics 
x <- bioticfit 
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- as.array(x[["Beta"]])[,,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(runs, 4, time)) # "4" is number of species

for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,,]
  n[z,,1] <- abundances[z,] #interations, species, time
  sigma <- diag(sigmas[z,])
  
  
  for(t in 2:time){
    
    #Biotic only
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1],
                             Sigma = sigma)
  }
}

#Create datafram for model checking
modelcheck_2022 <- n

#Create dataframe for plotting
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
abundances <- x[["n"]][,,14] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time)) #iterations, species, time


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,,]
  n[z,,1] <- abundances[z,]
  sigma <- diag(sigmas[z,])
  
  for(t in 2:time){
    
    #Remove env drivers
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1],
                             Sigma = sigma)
    
  }
}

#Create datafram for model checking
modelcheck_2023 <- n

sims2023median <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2023), c(2,3), median)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance") %>% 
  mutate(real_week = time + 24, year = 2023) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))
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
  #dplyr::mutate(real_week = time + 25, year = 2023) %>% 
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
n <- array(NA, dim = c(runs, 4, time)) #iterations, species, time


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,,]
  n[z,,1] <- abundances[z,]
  sigma <- diag(sigmas[z,])
  
  for(t in 2:time){
    
    #Remove env drivers
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1],
                             Sigma = sigma)
    
  }
}

#Create datafram for model checking
modelcheck_2024 <- n

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
simsbiotic <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(median = Abundance)

#Plot simulated predictions against latent states
# Plot 1: Only show Green Algae + Other N Fixers
p3 <- ggplot(simsbiotic, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsbiotic,
                               median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsbiotic,
                             median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), 
                                           median, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,59)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances: Only Biotic Interactions") +
  colScale + filScale + linScale + theme_bw()

# Plot 2: Only show Anabaena + Microcoleus
p4 <- ggplot(simsbiotic, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsbiotic,
                               median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsbiotic,
                             median = ifelse(Species %in% c("Anabaena", "Microcoleus"), 
                                           median, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 5)) +
  coord_cartesian(ylim = c(0,16)) +
  labs(x = "Date", y = "Percent Cover (%)") +
  colScale + filScale + linScale + theme_bw()


# Combine plots and collect legends
(p3 / p4) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

#Compile model check dataframes into a single full timeseries matrix
predictives.river.biotic <- abind(exp(modelcheck_2022), exp(modelcheck_2023), 
                                  exp(modelcheck_2024), along = 3)

#Save predictive output of Biotic model
saveRDS(predictives.river.biotic, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Biotic.rds"))


#-----------------------------------------------------------------------------
#River-wide: Abiotic Variables

#Pull out community abundances and demographics 
x <- abioticfit
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

Ptheta <- x[["Ptheta"]][,]
phos <- stand_nut$oPhos_ug_P_L[1:time]

Atheta <- x[["Atheta"]][,]
amon <- stand_nut$ammonium_mg_N_L[1:time]

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
  
  # #Pull env covariates
  nTheta <- Ntheta[z,]
  pTheta <- Ptheta[z,]
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      # Remove biotic interactions
        n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] + nTheta[s]*nitrate[t-1] +
                            pTheta[s]*phos[t-1] + aTheta[s]*amon[t-1] + 
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
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022median, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 25, year = 2022) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
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
Ntheta <- x[["Ntheta"]][,]
nitrate <- stand_nut$nitrate_mg_N_L[14:(13+time)]

Ptheta <- x[["Ptheta"]][,]
phos <- stand_nut$oPhos_ug_P_L[14:(13+time)]

Atheta <- x[["Atheta"]][,]
amon <- stand_nut$ammonium_mg_N_L[14:(13+time)]

Dtheta <- x[["Dtheta"]][,]
dis <- discharge$stand_discharge[14:(13+time)]

Ttheta <- x[["Ttheta"]][,]
temp <- stand_nut$temp_C[14:(13+time)]

Ctheta <- x[["Ctheta"]][,]
cond <- stand_nut$cond_uS_cm[14:(13+time)]

Rtheta <- x[["Rtheta"]][,]
rad <- swradiation$stand_rad[14:(13+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  # #Pull env covariates
  nTheta <- Ntheta[z,]
  pTheta <- Ptheta[z,]
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      # Remove biotic interactions
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] + nTheta[s]*nitrate[t-1] +
                          pTheta[s]*phos[t-1] + aTheta[s]*amon[t-1] + 
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
Ntheta <- x[["Ntheta"]][,]
nitrate <- stand_nut$nitrate_mg_N_L[29:(28+time)]

Ptheta <- x[["Ptheta"]][,]
phos <- stand_nut$oPhos_ug_P_L[29:(28+time)]

Atheta <- x[["Atheta"]][,]
amon <- stand_nut$ammonium_mg_N_L[29:(28+time)]

Dtheta <- x[["Dtheta"]][,]
dis <- discharge$stand_discharge[29:(28+time)]

Ttheta <- x[["Ttheta"]][,]
temp <- stand_nut$temp_C[29:(28+time)]

Ctheta <- x[["Ctheta"]][,]
cond <- stand_nut$cond_uS_cm[29:(28+time)]

Rtheta <- x[["Rtheta"]][,]
rad <- swradiation$stand_rad[29:(28+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,]
  n[z,,1] <- abundances[z,]
  sigma <- sigmas[z,]
  
  # #Pull env covariates
  nTheta <- Ntheta[z,]
  pTheta <- Ptheta[z,]
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      # Remove biotic interactions
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] + nTheta[s]*nitrate[t-1] +
                          pTheta[s]*phos[t-1] + aTheta[s]*amon[t-1] + 
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
simsabiotic <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(median = Abundance)

#Plot simulated predictions against latent states
# Plot 1: Only show Green Algae + Other N Fixers
p5 <- ggplot(simsabiotic, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsabiotic,
                               median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsabiotic,
                             median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), 
                                           median, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,75)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances: Only Abiotic Interactions") +
  colScale + filScale + linScale + theme_bw()

# Plot 2: Only show Anabaena + Microcoleus
p6 <- ggplot(simsabiotic, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsabiotic,
                               median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsabiotic,
                             median = ifelse(Species %in% c("Anabaena", "Microcoleus"), 
                                           median, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 5)) +
  coord_cartesian(ylim = c(0,16)) +
  labs(x = "Date", y = "Percent Cover (%)") +
  colScale + filScale + linScale + theme_bw()


# Combine plots and collect legends
(p5 / p6) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

#Compile model check dataframes into a single full timeseries matrix
predictives.river.abiotic <- abind(exp(modelcheck_2022), exp(modelcheck_2023), 
                                   exp(modelcheck_2024), along = 3)

#Save predictive output of Abiotic model
saveRDS(predictives.river.abiotic, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Abiotic.rds"))

#-----------------------------------------------------------------------------
#River-wide: Abiotic Variables Minus Nutrients

#Pull out community abundances and demographics 
x <- abioticnonutfit
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- x[["Beta"]][,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(runs, 4, time)) #iterations, species, time

#Pull out environmental effects
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
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    for(s in 1:4){

      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
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
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(exp(modelcheck_2022), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022median, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 25, year = 2022) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
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
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
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
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    for(s in 1:4){
      
      #Remove nutrients too
      n[z,s,t] <- rnorm(1, Alpha[s] + Beta[s]*n[z,s,t-1] +  
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
simsabioticnonut <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(median = Abundance)

#Plot simulated predictions against latent states
# Plot 1: Only show Green Algae + Other N Fixers
p7 <- ggplot(simsabioticnonut, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsabioticnonut,
                               median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsabioticnonut,
                             median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), 
                                           median, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,69)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances: Only Abiotic Interactions Minus Nutrients") +
  colScale + filScale + linScale + theme_bw()

# Plot 2: Only show Anabaena + Microcoleus
p8 <- ggplot(simsabioticnonut, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsabioticnonut,
                               median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsabioticnonut,
                             median = ifelse(Species %in% c("Anabaena", "Microcoleus"), 
                                           median, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 5)) +
  coord_cartesian(ylim = c(0,16)) +
  labs(x = "Date", y = "Percent Cover (%)") +
  colScale + filScale + linScale + theme_bw()


# Combine plots and collect legends
(p7 / p8) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

#Compile model check dataframes into a single full timeseries matrix
predictives.river.abioticnonut <- abind(exp(modelcheck_2022), exp(modelcheck_2023), 
                                   exp(modelcheck_2024), along = 3)

#Save predictive output of Abiotic Minus Nutrient model
saveRDS(predictives.river.abioticnonut, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_AbioticNoNut.rds"))

