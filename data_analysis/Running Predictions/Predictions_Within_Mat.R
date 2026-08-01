#Historical Predictions, 2022
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(abind)

#files mat_params2 and mat_params2_groups are in WithinMatModel_vs_real
#Read in cleaned latent states, params2_all. This file also reads in the cleaning_HAB.R file
source(here::here("data_analysis/Compare Obs Vs Modeled Outputs/Within_Mat_model_vs_real.R"))


#Graphing palettes
#Create a color palette
mycols <- c("brown", "#00538A", "#F6926A")
mypal <- palette(mycols)
names(mypal) <- c("Anabaena", "Epithemia Diatoms", "Geitlerinema")
colScale <- scale_color_manual(values = mypal)
filScale <- scale_fill_manual(values = mypal)
  
#Read in latent states and effect coefficients
M.fit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/WithinMat_Micro_predictions.rds"))

#Pull out community abundances and demographics 
x <- M.fit #m1 = averaged reaches for Microcoleus
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- as.array(x[["Beta"]])[,,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #13 weeks in 2022, 13 in 2023, 15 in 2024
n <- array(NA, dim = c(runs, 3, time)) #9 is number of species

#Pull out environmental effects
Ntheta <- x[["Ntheta"]][,]
nitrate <- stand_nut$nitrate_mg_N_L[1:time]

Ptheta <- x[["Ptheta"]][,]
phos <- stand_nut$oPhos_ug_P_L[1:time]

Atheta <- x[["Atheta"]][,]
amon <- stand_nut$ammonium_mg_N_L[1:time]

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

#Run predictions
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
  DINTheta <- DINtheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    #Include env drivers
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + 
                               DINTheta*DIN[t-1] +
                               dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)

  }
}

pred2022 <- n

sims2022median <- as.data.frame(t(as.data.frame(apply(exp(pred2022), c(2,3), median)))) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = c(1:3), names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(exp(pred2022), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(exp(pred2022), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIupper")

matsims2022 <- left_join(sims2022median, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))



#####################################
#Historical Predictions, 2023
#####################################
#Pull out community abundances and demographics 
abundances <- x[["n"]][,,14] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2023
n <- array(NA, dim = c(runs, 3, time))

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[14:(13+time)]

phos <- stand_nut$oPhos_ug_P_L[14:(13+time)]

amon <- stand_nut$ammonium_mg_N_L[14:(13+time)]

DIN <- stand_nut$DIN[14:(13+time)]

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
  DINTheta <- DINtheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  
  for(t in 2:time){
    #Everything included
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + 
                               DINTheta*DIN[t-1] +
                               dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    
    
  }
}

pred2023 <- n

sims2023median <- as.data.frame(t(as.data.frame(apply(exp(pred2023), c(2,3), median)))) %>% 
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = c(1:3), names_to = "Species", values_to = "Abundance")
sims2023lquant <- as.data.frame(t(as.data.frame(apply(exp(pred2023), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(exp(pred2023), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIupper")

matsims2023 <- left_join(sims2023median, sims2023lquant, by=c("Species", "time")) %>%
  left_join(., sims2023uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 26, year = 2023) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))




#####################################
#Historical Predictions, 2024
#####################################
#Pull out community abundances and demographics 
abundances <- x[["n"]][,,27] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(runs, 3, time))

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[27:(26+time)]

phos <- stand_nut$oPhos_ug_P_L[27:(26+time)]

amon <- stand_nut$ammonium_mg_N_L[27:(26+time)]

DIN <- stand_nut$DIN[27:(26+time)]

dis <- discharge$stand_discharge[27:(26+time)]

temp <- stand_nut$temp_C[27:(26+time)]

cond <- stand_nut$cond_uS_cm[27:(26+time)]

rad <- swradiation$stand_rad[27:(26+time)]


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
  DINTheta <- DINtheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    
    #Everything included
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + 
                               DINTheta*DIN[t-1] +
                               dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
  }
}

pred2024 <- n

sims2024median <- as.data.frame(t(as.data.frame(apply(exp(pred2024), c(2,3), median)))) %>% 
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = c(1:3), names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(exp(pred2024), c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(exp(pred2024), c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIupper")

matsims2024 <- left_join(sims2024median, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 26, year = 2024) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))


#Join together simulation data
matsimsallyears <- rbind(matsims2022, matsims2023, matsims2024) %>% 
  dplyr::rename(median = Abundance)

###Create plot of TM microscopy predictions vs latent states
ggplot(matsimsallyears, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species), alpha = 0.3) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5) +
  # Latent points/lines
  geom_line(data = mat_params2_TM, aes(linetype = "Latent", colour = Species), 
            linewidth = 2) +
  scale_y_continuous(breaks = seq(0, 600, 5)) +
  coord_cartesian(ylim = c(0,15)) +
  labs(x = "Date", y = "Relative Abundance (%)", title = "Within-Mat Latent vs. Predicted Abundances") +
  colScale + filScale + linScale + theme_bw()



#Compile model check dataframes into a single full timeseries matrix
predictives_TMmats <- abind(exp(pred2022), exp(pred2023), exp(pred2024), along = 3)

#Save predictive output of Microcoleus model
saveRDS(predictives_TMmats, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/WithinMat_Pred_TM.rds"))


