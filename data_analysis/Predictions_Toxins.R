#Historical Predictions of Anatoxin Concentrations, 2022
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(abind)

#files mat_params2 and mat_params2_groups are in WithinMatModel_vs_real

#Graphing palettes
#Create a color palette
mycols <- c("brown", "darkolivegreen4", "darkorange", "chartreuse3", 
            "lavender", "darkcyan", "mediumpurple3","khaki1", "antiquewhite3",
            "goldenrod", "lightblue1")
mypal <- palette(mycols)
names(mypal) = c("Anabaena", "Epithemia Diatoms", "Geitlerinema", 
                 "Green Algae", "Leptolyngbya", "Microcoleus", 
                 "Non-Epithemia Diatoms", "Nostoc", "Oscillatoria",
                 "Other Coccoids", "Rare")
colScale <- scale_color_manual(values = mypal)
filScale <- scale_fill_manual(values = mypal)

#Read in latent states and effect coefficients
M.fit <- readRDS(here::here("data/WithinMat_Micro_predictions.rds"))
A.fit <- readRDS(here::here("data/WithinMat_Ana_predictions.rds"))


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
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    #Include env drivers
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    
  }
}

pred2022 <- n

# 'Green Algae' = V4, 
# Microcoleus = V5, 'Non-Epithemia Diatoms' = V6,
# Nostoc = V7, 'Other Coccoids' = V8,
# Rare = V9)

sims2022mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:3, exp)) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = c(1:3), names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:3, exp)) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:3, exp)) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIupper")

matsims2022 <- left_join(sims2022mean, sims2022lquant, by=c("Species", "time")) %>%
  left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))


#Plot

mat_p22 <- ggplot(subset(matsims2023, Species %in% c("Anabaena",
                                                     "Epithemia Diatoms",
                                                     "Geitlerinema")),
                  aes(x = model_date, y = Abundance)) +
  geom_line(size = 1.5, aes(color = Species)) +
  geom_line(data = subset(mat_params2_TM[mat_params2_TM$year %in% "2023", ],
                          Species %in% c("Anabaena",
                                         "Epithemia Diatoms",
                                         "Geitlerinema")),
            aes(x = model_date, y = mean, colour = Species),
            linewidth = 4, alpha = .35) +
  scale_y_continuous(breaks=c(seq(0,150,5))) +
  #coord_cartesian(ylim = c(0,70)) +
  labs(x = "Date", y = "Relative Abundance (%)", title = "2023 Predictions") +
  colScale


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
    
    
  }
}

pred2023 <- n

sims2023mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:3, exp)) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = c(1:3), names_to = "Species", values_to = "Abundance")
sims2023lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:3, exp)) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIlower")
sims2023uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:3, exp)) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:3, names_to = "Species", values_to = "CIupper")

matsims2023 <- left_join(sims2023mean, sims2023lquant, by=c("Species", "time")) %>%
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
  }
}

pred2024 <- n

sims2024mean <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:9, exp)) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3, 'Green Algae' = V4, 
                Microcoleus = V5, 'Non-Epithemia Diatoms' = V6,
                Nostoc = V7, 'Other Coccoids' = V8,
                Rare = V9) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = c(1:9), names_to = "Species", values_to = "Abundance")
sims2024lquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.025)))) %>% 
  dplyr::mutate(across(1:9, exp)) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3, 'Green Algae' = V4, 
                Microcoleus = V5, 'Non-Epithemia Diatoms' = V6,
                Nostoc = V7, 'Other Coccoids' = V8,
                Rare = V9) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:9, names_to = "Species", values_to = "CIlower")
sims2024uquant <- as.data.frame(t(as.data.frame(apply(n, c(2,3), quantile, probs = 0.975)))) %>% 
  dplyr::mutate(across(1:9, exp)) %>%
  dplyr::rename(Anabaena = V1, 'Epithemia Diatoms' = V2,
                Geitlerinema = V3, 'Green Algae' = V4, 
                Microcoleus = V5, 'Non-Epithemia Diatoms' = V6,
                Nostoc = V7, 'Other Coccoids' = V8,
                Rare = V9) %>% 
  dplyr::mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:9, names_to = "Species", values_to = "CIupper")

matsims2024 <- left_join(sims2024mean, sims2024lquant, by=c("Species", "time")) %>%
  left_join(., sims2024uquant, by=c("Species", "time")) %>% 
  dplyr::mutate(real_week = time + 26, year = 2024) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))


#Join together simulation data
matsimsallyears <- rbind(matsims2022, matsims2023, matsims2024) %>% 
  dplyr::rename(mean = Abundance)

ggplot(subset(matsims2022, Species %in% c("Anabaena",
                                          "Epithemia Diatoms",
                                          "Geitlerinema")),
       aes(x = model_date, y = Abundance)) +
  geom_line(size = 1.5, aes(color = Species)) +
  geom_line(data = subset(mat_params2_TM[mat_params2_TM$year %in% "2022", ],
                          Species %in% c("Anabaena",
                                         "Epithemia Diatoms",
                                         "Geitlerinema")),
            aes(x = model_date, y = mean, colour = Species),
            linewidth = 4, alpha = .35) +
  scale_y_continuous(breaks=c(seq(0,150,5))) +
  #coord_cartesian(ylim = c(0,70)) +
  labs(x = "Date", y = "Relative Abundance (%)", title = "2022 Predictions") +
  colScale

###Create plot of TM microscopy predictions vs latent states
ggplot(subset(matsimsallyears, Species %in% c("Anabaena", "Epithemia Diatoms",
                                              "Geitlerinema")), 
       aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species), alpha = 0.3) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", colour = Species), size = 1.5) +
  # Latent points/lines
  geom_line(data = subset(mat_params2_TM, Species %in% c("Anabaena", "Epithemia Diatoms",
                                                         "Geitlerinema")),
            aes(linetype = "Latent", colour = Species), linewidth = 2) +
  scale_y_continuous(breaks = seq(0, 600, 10)) +
  coord_cartesian(ylim = c(0,20)) +
  labs(x = "Date", y = "Relative Abundance (%)", title = "Latent vs. Predicted Abundances") +
  colScale + filScale + linScale



#Compile model check dataframes into a single full timeseries matrix
predictives_TMmats <- abind(pred2022, pred2023, pred2024, along = 3)
predictives_TMmats <- exp(predictives_TMmats)

#Save predictive output of Microcoleus model
saveRDS(predictives_TMmats, 
        file = here::here("data/WithinMat_Pred_TM.rds"))