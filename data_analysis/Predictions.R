#Historical Predictions, 2022
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(abind)

#Read in latent states, params2_all
#source(here::here("data_analysis/model_vs_real_data.R"))

#Create plot palettes for graphing
#Create a color palette
mycols <- c("brown", "darkolivegreen4", "darkcyan", "darkorange")
mypal <- palette(mycols)
names(mypal) = c("Anabaena", "Green Algae", "Microcoleus", 
                 "Other N Fixers")
colScale <- scale_color_manual(values = mypal)
filScale <- scale_fill_manual(values = mypal)
linScale <- scale_linetype_manual(name = "State Type",
                                  values = c("Latent" = "11",
                                             "Predicted" = "solid"))

#Read in latent states and effect coefficients
allfit <- readRDS(here::here("data/Riverwide_AllVar_predictions.rds"))
bioticfit <- readRDS(here::here("data/Riverwide_Biotic_predictions.rds"))
abioticfit <- readRDS(here::here("data/Riverwide_Abiotic_predictions.rds"))
abioticnonutfit <- readRDS(here::here("data/Riverwide_AbioticNonut_predictions.rds"))

#----------------------------------------------------------------------------
#River-wide: All Variables

#Pull out community abundances and demographics 
x <- allfit #m4 = all vars, m5 = biotic only
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- as.array(x[["Beta"]])[,,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(runs, 4, time)) # "4" is number of species

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
  Beta <- betas[z,,]
  n[z,,1] <- abundances[z,] #interations, species, time
  sigma <- diag(sigmas[z,])
  
  #Pull env covariates
  nTheta <- Ntheta[z,]
  pTheta <- Ptheta[z,]
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
      
    #Everything included
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    # #Biotic only
    # n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1],
    #                          Sigma = sigma)
  }
}

#Create datafram for model checking
modelcheck_2022 <- n

#Create dataframe for plotting
sims2022mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022mean, sims2022lquant, by=c("Species", "time")) %>%
                        left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))

# #Plot
# ggplot(sims2022, aes(x = model_date, y = Abundance, colour = Species)) +
#   geom_line(aes(linetype = "Predicted"), size = 1.5) +
#   geom_line(data=params2_all[params2_all$year %in% "2022", ],
#             aes(x = model_date, y = mean, colour = Species, linetype = "Latent"),
#             linewidth = 2) +
#   geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`,
#                   fill = Species), color = NA, alpha = 0.2) +
#   scale_y_continuous(breaks=c(seq(0,810,30))) +
#   coord_cartesian(ylim = c(0,810)) +
#   labs(x = "Date", y = "Percent Cover (%)", title = "2022 Predictions") +
#   colScale + filScale + linScale
  



#
#Historical Predictions, 2023
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,14] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time)) #iterations, species, time

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[14:(13+time)]

phos <- stand_nut$oPhos_ug_P_L[14:(13+time)]

amon <- stand_nut$ammonium_mg_N_L[14:(13+time)]

dis <- discharge$stand_discharge[14:(13+time)]

temp <- stand_nut$temp_C[14:(13+time)]

cond <- stand_nut$cond_uS_cm[14:(13+time)]

rad <- swradiation$stand_rad[14:(13+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,,]
  n[z,,1] <- abundances[z,]
  sigma <- diag(sigmas[z,])
  
  #Pull env covariates
  nTheta <- Ntheta[z,]
  pTheta <- Ptheta[z,]
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    
    #Everything included
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    
    # #Remove env drivers
    # n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1],
    #                          Sigma = sigma)

    
  }
}

#Create datafram for model checking
modelcheck_2023 <- n

sims2023mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance") %>% 
  mutate(real_week = time + 24, year = 2023) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))
sims2023lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2023 <- left_join(sims2023mean, sims2023lquant, by=c("Species", "time")) %>%
  left_join(., sims2023uquant, by=c("Species", "time")) %>% 
  #dplyr::mutate(real_week = time + 25, year = 2023) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))



# p23 <- ggplot(sims2023, aes(x = model_date, y = Abundance, colour = Species)) +
#   geom_line(aes(linetype = "Predicted"), size = 1.5) +
#   geom_line(data=params2_all[params2_all$year %in% "2023", ], 
#             aes(x = model_date, y = mean, colour = Species, linetype = "Latent"), 
#             linewidth = 2) +
#   geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`,
#                   fill = Species), color = NA, alpha = 0.2) +
#   scale_y_continuous(breaks=c(seq(0,100,5))) +
#   coord_cartesian(ylim = c(0,50)) +
#   labs(x = "Date", y = "Percent Cover (%)", title = "2023 Predictions") +
#   colScale + filScale + linScale


#
#Historical Predictions, 2024
#

#Pull out community abundances and demographics 
abundances <- x[["n"]][,,29] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 17 #number of weeks in 2023
n <- array(NA, dim = c(runs, 4, time)) #iterations, species, time

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[29:(28+time)]

phos <- stand_nut$oPhos_ug_P_L[29:(28+time)]

amon <- stand_nut$ammonium_mg_N_L[29:(28+time)]

dis <- discharge$stand_discharge[29:(28+time)]

temp <- stand_nut$temp_C[29:(28+time)]

cond <- stand_nut$cond_uS_cm[29:(28+time)]

rad <- swradiation$stand_rad[29:(28+time)]


for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,,]
  n[z,,1] <- abundances[z,]
  sigma <- diag(sigmas[z,])
  
  #Pull env covariates
  nTheta <- Ntheta[z,]
  pTheta <- Ptheta[z,]
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    
    #Everything included
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    # #Remove env drivers
    # n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1],
    #                          Sigma = sigma)
    
  }
}

#Create datafram for model checking
modelcheck_2024 <- n

sims2024mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2024 <- left_join(sims2024mean, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 24, year = 2024) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))


# p24 <- ggplot(sims2024, aes(x = model_date, y = Abundance, colour = Species)) +
#   geom_line(aes(linetype = "Predicted"), size = 1.5) +
#   geom_line(data=params2_all[params2_all$year %in% "2024", ],
#             aes(x = model_date, y = mean, colour = Species, linetype = "Latent"),
#             linewidth = 2) +
#   geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`,
#                   fill = Species), color = NA, alpha = 0.2) +
#   scale_y_continuous(breaks=c(seq(0,100,5))) +
#   coord_cartesian(ylim = c(0,50)) +
#   labs(x = "Date", y = "Percent Cover (%)", title = "2023 Predictions") +
#   colScale + filScale + linScale
# 
# 
# ggarrange(
#   p22, p23, p24, labels = c("A", "B", "C"), ncol = 3,
#   common.legend = TRUE, legend = "bottom"
# )

#Join together simulation data
simsallyears <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(mean = Abundance)

#Plot simulated predictions against latent states
# Plot 1: Only show Green Algae + Other N Fixers
p1 <- ggplot(simsallyears, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsallyears,
                               mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), mean, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
             data = transform(simsallyears,
                              mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), 
                                            mean, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), mean, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,100)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances") +
  colScale + filScale + linScale

# Plot 2: Only show Anabaena + Microcoleus
p2 <- ggplot(simsallyears, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsallyears,
                               mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), mean, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsallyears,
                             mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), 
                                           mean, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), mean, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,60)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances") +
  colScale + filScale + linScale
  

# Combine plots and collect legends
(p1 / p2) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical")


#Compile model check dataframes into a single full timeseries matrix
predictives.river.all <- abind(modelcheck_2022, modelcheck_2023, 
                        modelcheck_2024, along = 3)
predictives.river.all <- exp(predictives.river.all)

#Save predictive output of All Variables model
saveRDS(predictives.river.all, 
        file = here::here("data/Riverwide_Pred_AllVar.rds"))



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
sims2022mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022mean, sims2022lquant, by=c("Species", "time")) %>%
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

sims2023mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance") %>% 
  mutate(real_week = time + 24, year = 2023) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))
sims2023lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2023 <- left_join(sims2023mean, sims2023lquant, by=c("Species", "time")) %>%
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

sims2024mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2024 <- left_join(sims2024mean, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 24, year = 2024) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))

#Join together simulation data
simsallyears <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(mean = Abundance)

#Plot simulated predictions against latent states
# Plot 1: Only show Green Algae + Other N Fixers
p3 <- ggplot(simsallyears, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsallyears,
                               mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), mean, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsallyears,
                             mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), 
                                           mean, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), mean, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,70)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances: Only Biotic Interactions") +
  colScale + filScale + linScale

# Plot 2: Only show Anabaena + Microcoleus
p4 <- ggplot(simsallyears, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsallyears,
                               mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), mean, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsallyears,
                             mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), 
                                           mean, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), mean, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,20)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances") +
  colScale + filScale + linScale


# Combine plots and collect legends
(p3 / p4) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

#Compile model check dataframes into a single full timeseries matrix
predictives.river.biotic <- abind(modelcheck_2022, modelcheck_2023, 
                               modelcheck_2024, along = 3)
predictives.river.biotic <- exp(predictives.river.biotic)

#Save predictive output of All Variables model
saveRDS(predictives.river.biotic, 
        file = here::here("data/Riverwide_Pred_Biotic.rds"))


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
        n[z,s,t] <- Alpha[s] + Beta[s]*n[z,s,t-1] + nTheta[s]*nitrate[t-1]+
          pTheta[s]*phos[t-1] + aTheta[s]*amon[t-1] + dTheta[s]*dis[t-1] +
          tTheta[s]*temp[t-1] + cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1]
    }
  }
}

#Create datafram for model checking
modelcheck_2022 <- n

#Create dataframe
sims2022mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022mean, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 25, year = 2022) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



#Plot
ggplot(sims2022, aes(x = model_date, y = Abundance, colour = Species)) +
  geom_line(aes(linetype = "Predicted"), size = 1.5) +
  geom_line(data=params2_all[params2_all$year %in% "2022", ],
            aes(x = model_date, y = mean, colour = Species, linetype = "Latent"),
            linewidth = 2) +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`,
                  fill = Species), color = NA, alpha = 0.2) +
  scale_y_continuous(breaks=c(seq(0,810,30))) +
  coord_cartesian(ylim = c(0,810)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "2022 Predictions") +
  colScale + filScale + linScale
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
      n[z,s,t] <- Alpha[s] + Beta[s]*n[z,s,t-1] + nTheta[s]*nitrate[t-1]+
        pTheta[s]*phos[t-1] + aTheta[s]*amon[t-1] + dTheta[s]*dis[t-1] +
        tTheta[s]*temp[t-1] + cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1]
    }
  }
}

#Create datafram for model checking
modelcheck_2023 <- n

#Create dataframe
sims2023mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2023lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2023 <- left_join(sims2023mean, sims2023lquant, by=c("Species", "time")) %>%
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
      n[z,s,t] <- Alpha[s] + Beta[s]*n[z,s,t-1] + nTheta[s]*nitrate[t-1]+
        pTheta[s]*phos[t-1] + aTheta[s]*amon[t-1] + dTheta[s]*dis[t-1] +
        tTheta[s]*temp[t-1] + cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1]
    }
  }
}

#Create datafram for model checking
modelcheck_2024 <- n

#Create dataframe
sims2024mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2024 <- left_join(sims2024mean, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 24, year = 2024) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



#Join together simulation data
simsallyears <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(mean = Abundance)

#Plot simulated predictions against latent states
# Plot 1: Only show Green Algae + Other N Fixers
p5 <- ggplot(simsallyears, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsallyears,
                               mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), mean, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsallyears,
                             mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), 
                                           mean, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), mean, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,80)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances: Only Abiotic Interactions Minus Nutrients") +
  colScale + filScale + linScale

# Plot 2: Only show Anabaena + Microcoleus
p6 <- ggplot(simsallyears, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsallyears,
                               mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), mean, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsallyears,
                             mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), 
                                           mean, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), mean, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,20)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances") +
  colScale + filScale + linScale


# Combine plots and collect legends
(p5 / p6) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

#Compile model check dataframes into a single full timeseries matrix
predictives.river.abiotic <- abind(modelcheck_2022, modelcheck_2023, 
                                        modelcheck_2024, along = 3)
predictives.river.abiotic <- exp(predictives.river.abiotic)

#Save predictive output of All Variables model
saveRDS(predictives.river.abiotic, 
        file = here::here("data/Riverwide_Pred_Abiotic.rds"))

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
        n[z,s,t] <- Alpha[s] + Beta[s]*n[z,s,t-1] + dTheta[s]*dis[t-1] +
          tTheta[s]*temp[t-1] + cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1]
    }
  }
}

#Create datafram for model checking
modelcheck_2022 <- n

#Create dataframe
sims2022mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022mean, sims2022lquant, by=c("Species", "time")) %>%
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
      n[z,s,t] <- Alpha[s] + Beta[s]*n[z,s,t-1] + dTheta[s]*dis[t-1] +
        tTheta[s]*temp[t-1] + cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1]
      
    }
  }
}

#Create datafram for model checking
modelcheck_2023 <- n

#Create dataframe
sims2023mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2023lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2023 <- left_join(sims2023mean, sims2023lquant, by=c("Species", "time")) %>%
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
      n[z,s,t] <- Alpha[s] + Beta[s]*n[z,s,t-1] + dTheta[s]*dis[t-1] +
        tTheta[s]*temp[t-1] + cTheta[s]*cond[t-1] + rTheta[s]*rad[t-1]
      
    }
  }
}

#Create datafram for model checking
modelcheck_2024 <- n

#Create dataframe
sims2024mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:4, exp)) %>%
  dplyr::rename("Green Algae" = V1, "Microcoleus" = V2,
                "Anabaena" = V3,
                "Other N Fixers" = V4) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2024 <- left_join(sims2024mean, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 24, year = 2024) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



#Join together simulation data
simsallyears <- rbind(sims2022, sims2023, sims2024) %>% 
  dplyr::rename(mean = Abundance)

#Plot simulated predictions against latent states
# Plot 1: Only show Green Algae + Other N Fixers
p7 <- ggplot(simsallyears, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsallyears,
                               mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), mean, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsallyears,
                             mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), 
                                           mean, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), mean, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,80)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances: Only Abiotic Interactions Minus Nutrients") +
  colScale + filScale + linScale

# Plot 2: Only show Anabaena + Microcoleus
p8 <- ggplot(simsallyears, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(simsallyears,
                               mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), mean, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5,
            data = transform(simsallyears,
                             mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), 
                                           mean, NA))) +
  # Latent points/lines
  geom_line(aes(linetype = "Latent", colour = Species), linewidth = 2,
            data = transform(params2_all, mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), mean, NA))) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,20)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Latent vs. Predicted Abundances") +
  colScale + filScale + linScale


# Combine plots and collect legends
(p7 / p8) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical")

#Compile model check dataframes into a single full timeseries matrix
predictives.river.abioticnonut <- abind(modelcheck_2022, modelcheck_2023, 
                                   modelcheck_2024, along = 3)
predictives.river.abioticnonut <- exp(predictives.river.abioticnonut)

#Save predictive output of All Variables model
saveRDS(predictives.river.abioticnonut, 
        file = here::here("data/Riverwide_Pred_AbioticNoNut.rds"))

