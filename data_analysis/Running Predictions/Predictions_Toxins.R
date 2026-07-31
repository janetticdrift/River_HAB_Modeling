#Historical Predictions of Anatoxin Concentrations, 2022
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(abind)

#Read in raw, uncleaned latent states and effect coefficients for toxin concentrations (River.fit and
  #Mat.fit),
  #and latent states of percent cover and microscopy abundances 
  #(percentcover_latent and microscopyTM_latent)
River.fit.TM <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_predictions.rds"))
River.fit.TAC <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TAC_River_predictions.rds"))
Mat.fit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_Mat_predictions.rds"))
percentcover_latent <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_LatentStates.rds"))
microscopyTM_latent <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/WithinMat_Micro_LatentStates.rds"))

#Read in cleaned toxin latent states, used for creating graphs
source(here::here("data_analysis/Compare Obs Vs Modeled Outputs/Toxins_model_vs_real.R"))

######Simulate Toxins using microscopy data------

#Pull out community abundances and demographics 
x <- Mat.fit 
tox_conc <- x[["tox_raw"]][,1:2] #iterations, time
Beta0 <- x[["Beta0"]]
Beta1 <- x[["Beta1"]]
Beta2 <- x[["Beta2"]]
Beta3 <- x[["Beta3"]]
sigma_p <- x[["sigma_p"]]

#Pull out environmental effects
Ntheta <- x[["Ntheta"]]
Ptheta <- x[["Ptheta"]]
Atheta <- x[["Atheta"]]
#DINtheta <- x[["DINtheta"]]
Dtheta <- x[["Dtheta"]]
Ttheta <- x[["Ttheta"]]
Ctheta <- x[["Ctheta"]]
Rtheta <- x[["Rtheta"]]

#Create design matrix
X1 <- cbind(
  intercept = rep(1, time),
  microscopyTM_latent[, -1][1:time, ],  #Abundances are on the log-scale
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][1:time],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][1:time],
  amon = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][1:time],
  #DIN = stand_nut$DIN[-c(14:15, 29:30)][1:time],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)][1:time],
  temp = stand_nut$temp_C[-c(14:15, 29:30)][1:time],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)][1:time],
  rad = swradiation$stand_rad[-c(14:15, 29:30)][1:time]
)

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 13

tox <- matrix(NA, runs, time)

for (z in 1:runs) {
  
  # Build parameter vectors
  beta <- c(
    Beta0[z],
    Beta1[z],
    Beta2[z],
    Beta3[z],
    Ntheta[z],
    Ptheta[z],
    Atheta[z],
    #DINtheta[z],
    Dtheta[z],
    Ttheta[z],
    Ctheta[z],
    Rtheta[z]
  )
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  
  #Simulation
  for (t in 2:time) {
    
    tox[z,t] <- rnorm(1, sum(X1[t-1, ] * beta), sigma_p[z])
  }
}

pred2022 <- tox

sims2022median <- as.data.frame(apply(exp(pred2022)/1000, 2, median)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(time = 1:time) 
sims2022lquant <- as.data.frame(apply(exp(pred2022)/1000, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>%
  dplyr::mutate(time = 1:time)
sims2022uquant <- as.data.frame(apply(exp(pred2022)/1000, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(time = 1:time)

#Join together median and CI ranges
matsims2022 <- dplyr::left_join(sims2022median, sims2022lquant, by=c("time")) %>%
  dplyr::left_join(., sims2022uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))


####
#Historical Predictions, 2023
####
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,14:15] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 13
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- cbind(
  intercept = rep(1, time),
  microscopyTM_latent[, -1][14:(13+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][14:(13+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][14:(13+time)],
  amon = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][14:(13+time)],
  #DIN = stand_nut$DIN[-c(14:15, 29:30)][14:(13+time)],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)][14:(13+time)],
  temp = stand_nut$temp_C[-c(14:15, 29:30)][14:(13+time)],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)][14:(13+time)],
  rad = swradiation$stand_rad[-c(14:15, 29:30)][14:(13+time)]
)


for (z in 1:runs) {
  
  # Build parameter vectors
  beta <- c(
    Beta0[z],
    Beta1[z],
    Beta2[z],
    Beta3[z],
    Ntheta[z],
    Ptheta[z],
    Atheta[z],
    #DINtheta[z],
    Dtheta[z],
    Ttheta[z],
    Ctheta[z],
    Rtheta[z]
  )
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  
  #Simulation
  for (t in 2:time) {
    
    tox[z,t] <- rnorm(1, sum(X1[t-1, ] * beta), sigma_p[z])
  }
}

pred2023 <- tox

sims2023median <- as.data.frame(apply(exp(pred2023)/1000, 2, median)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(time = 1:time) 
sims2023lquant <- as.data.frame(apply(exp(pred2023)/1000, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(time = 1:time)
sims2023uquant <- as.data.frame(apply(exp(pred2023)/1000, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(time = 1:time)

#Join together median and CI ranges
matsims2023 <- dplyr::left_join(sims2023median, sims2023lquant, by=c("time")) %>%
  dplyr::left_join(., sims2023uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 26, year = 2023) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))



###
#Historical Predictions, 2024
###
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,27:28] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 15
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- cbind(
  intercept = rep(1, time),
  microscopyTM_latent[, -1][27:(26+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][27:(26+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][27:(26+time)],
  amon = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][27:(26+time)],
  #DIN = stand_nut$DIN[-c(14:15, 29:30)][27:(26+time)],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)][27:(26+time)],
  temp = stand_nut$temp_C[-c(14:15, 29:30)][27:(26+time)],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)][27:(26+time)],
  rad = swradiation$stand_rad[-c(14:15, 29:30)][27:(26+time)]
)

for (z in 1:runs) {
  
  # Build parameter vectors
  beta <- c(
    Beta0[z],
    Beta1[z],
    Beta2[z],
    Beta3[z],
    Ntheta[z],
    Ptheta[z],
    Atheta[z],
    #DINtheta[z],
    Dtheta[z],
    Ttheta[z],
    Ctheta[z],
    Rtheta[z]
  )
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  
  #Simulation
  for (t in 2:time) {
    
    tox[z,t] <- rnorm(1, sum(X1[t-1, ]*beta), sigma_p[z])
  }
}

pred2024 <- tox

sims2024median <- as.data.frame(apply(exp(pred2024)/1000, 2, median)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(time = 1:time) 
sims2024lquant <- as.data.frame(apply(exp(pred2024)/1000, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(time = 1:time)
sims2024uquant <- as.data.frame(apply(exp(pred2024)/1000, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(time = 1:time)

#Join together median and CI ranges
matsims2024 <- dplyr::left_join(sims2024median, sims2024lquant, by=c("time")) %>%
  dplyr::left_join(., sims2024uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 26, year = 2024) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))

#Join together simulation data
matsimsallyears <- rbind(matsims2022, matsims2023, matsims2024)

###Create plot of TM microscopy predictions vs latent states

#Graphing palettes
#Create a color palette
mycols <- c("#791c55", "#af6c5d")
mypal <- palette(mycols)
mypal <- palette(mycols)
names(mypal) = c("Latent", "Predicted")
matcolScale <- scale_color_manual(name = "State Type", values = mypal)
filScale <- scale_fill_manual(name = "State Type", values = mypal)
linScale <- scale_linetype_manual(name = "State Type",
                                  values = c("Latent" = "11",
                                             "Predicted" = "solid"))

ggplot(matsimsallyears, aes(x = model_date, y = toxins)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = "Predicted"), 
              alpha = 0.3) +
  # Latent points/lines
  geom_line(data = tox_params2_mat,
            aes(y = median, linetype = "Latent", color = "Latent"), linewidth = 2) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", color = "Predicted"), size = 1.5) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(ylim = c(0,40)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "Within-Mat: Latent vs. Predicted Toxin Concentrations") +
  matcolScale + linScale + filScale + theme_bw()



#Compile model check dataframes into a single full timeseries matrix
predictives_toxins_mats1 <- abind(pred2022, pred2023, pred2024, along = 2)
predictives_toxins_mats <- exp(predictives_toxins_mats1)

#Save predictive output of Toxin Within-Mat model
saveRDS(predictives_toxins_mats, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Toxins_Pred_Withinmat.rds")) 




######Simulate Microcoleus Toxins using percent cover data-----

#Pull out community abundances and demographics 
x <- River.fit.TM 
tox_conc <- x[["tox_raw"]][,1:2] #iterations, time
Beta0 <- x[["Beta0"]]
Beta1 <- x[["Beta1"]]
Beta2 <- x[["Beta2"]]
Beta3 <- x[["Beta3"]]
Beta4 <- x[["Beta4"]]
BetaGreen <- x[["BetaGreen"]]
BetaMicro <- x[["BetaMicro"]]
BetaAna <- x[["BetaAna"]]
BetaNFix <- x[["BetaNFix"]]

sigma_p <- x[["sigma_p"]]

#Pull out environmental effects
Ntheta <- x[["Ntheta"]]
Ptheta <- x[["Ptheta"]]
Atheta <- x[["Atheta"]]
Dtheta <- x[["Dtheta"]]
Ttheta <- x[["Ttheta"]]
Ctheta <- x[["Ctheta"]]
Rtheta <- x[["Rtheta"]]

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 13

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  percentcover_latent[-c(14:15, 29:30), -1][1:time, ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][1:time],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][1:time],
  amon = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][1:time],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)][1:time],
  temp = stand_nut$temp_C[-c(14:15, 29:30)][1:time],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)][1:time],
  rad = swradiation$stand_rad[-c(14:15, 29:30)][1:time]
))

tox <- matrix(NA, runs, time)

# Build parameter matrixes
beta_matrix <- cbind(Beta0, Beta1, Beta2, Beta3, Beta4, Ntheta, Ptheta, Atheta, Dtheta, Ttheta,
              Ctheta, Rtheta)
beta_lag_matrix <- cbind(BetaGreen, BetaMicro, BetaAna, BetaNFix)

for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  beta_lag <- beta_lag_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, beta%*%X1[t-1, ] + beta_lag%*%X1[t-2, 2:5],
                         sigma_p[z])
  }
}

pred2022 <- tox

sims2022median <- as.data.frame(apply(exp(pred2022)/1000, 2, median)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(time = 1:time) 
sims2022lquant <- as.data.frame(apply(exp(pred2022)/1000, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(time = 1:time)
sims2022uquant <- as.data.frame(apply(exp(pred2022)/1000, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(time = 1:time)

#Join together median and CI ranges
riversims2022 <- dplyr::left_join(sims2022median, sims2022lquant, by=c("time")) %>%
  dplyr::left_join(., sims2022uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7))


###
#Historical Predictions, 2023
###
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,14:15] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 13
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  percentcover_latent[-c(14:15, 29:30), -1][14:(13+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][14:(13+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][14:(13+time)],
  amon = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][14:(13+time)],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)][14:(13+time)],
  temp = stand_nut$temp_C[-c(14:15, 29:30)][14:(13+time)],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)][14:(13+time)],
  rad = swradiation$stand_rad[-c(14:15, 29:30)][14:(13+time)]
))


for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  beta_lag <- beta_lag_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, beta%*%X1[t-1, ] + beta_lag%*%X1[t-2, 2:5],
                         sigma_p[z])
  }
}

pred2023 <- tox

sims2023median <- as.data.frame(apply(exp(pred2023)/1000, 2, median)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(time = 1:time) 
sims2023lquant <- as.data.frame(apply(exp(pred2023)/1000, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(time = 1:time)
sims2023uquant <- as.data.frame(apply(exp(pred2023)/1000, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(time = 1:time)

#Join together median and CI ranges
riversims2023 <- dplyr::left_join(sims2023median, sims2023lquant, by=c("time")) %>%
  dplyr::left_join(., sims2023uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 26, year = 2023) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))



###
#Historical Predictions, 2024
###
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,27:28] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 15
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  percentcover_latent[, -1][27:(26+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][27:(26+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][27:(26+time)],
  amon = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][27:(26+time)],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)][27:(26+time)],
  temp = stand_nut$temp_C[-c(14:15, 29:30)][27:(26+time)],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)][27:(26+time)],
  rad = swradiation$stand_rad[-c(14:15, 29:30)][27:(26+time)]
))

for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  beta_lag <- beta_lag_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, beta%*%X1[t-1, ] + beta_lag%*%X1[t-2, 2:5],
                         sigma_p[z])
  }
}

pred2024 <- tox

sims2024median <- as.data.frame(apply(exp(pred2024)/1000, 2, median)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(time = 1:time) 
sims2024lquant <- as.data.frame(apply(exp(pred2024)/1000, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(time = 1:time)
sims2024uquant <- as.data.frame(apply(exp(pred2024)/1000, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(time = 1:time)

#Join together median and CI ranges
riversims2024 <- dplyr::left_join(sims2024median, sims2024lquant, by=c("time")) %>%
  dplyr::left_join(., sims2024uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 26, year = 2024) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))

#Join together simulation data
riverTMsimsallyears <- rbind(riversims2022, riversims2023, riversims2024)

###Create plot of River-wide predictions vs latent states

#Graphing palettes
#Create a color palette
mycols <- c("#791C55", "#41789A")
mypal <- palette(mycols)
mypal <- palette(mycols)
names(mypal) = c("Latent", "Predicted")
rivercolScale <- scale_color_manual(name = "State Type", values = mypal)
filScale <- scale_fill_manual(name = "State Type", values = mypal)
linScale <- scale_linetype_manual(name = "State Type",
                                  values = c("Latent" = "11",
                                             "Predicted" = "solid"))

ggplot(riverTMsimsallyears, aes(x = model_date, y = toxins)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = "Predicted"), 
              alpha = 0.3) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", color = "Predicted"), size = 1.5) +
  # Latent points/lines
  geom_line(data = tox_params2_riverTM,
            aes(y = median, linetype = "Latent", color = "Latent"), linewidth = 2) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(ylim = c(0,40)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "River-Wide Toxin") +
  rivercolScale + filScale + linScale + theme_bw()



#Compile model check dataframes into a single full timeseries matrix
predictives_TM_toxins_river <- abind(pred2022, pred2023, pred2024, along = 2)
predictives_TM_toxins_river <- exp(predictives_TM_toxins_river)

#Save predictive output of Toxin River-Wide model
saveRDS(predictives_TM_toxins_river, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Toxins_TM_Pred_RiverWide.rds")) 




######Simulate Anabaena Toxins using percent cover data-----

#Pull out community abundances and demographics 
x <- River.fit.TAC 
tox_conc <- x[["tox_raw"]][,1:2] #iterations, time
Beta0 <- x[["Beta0"]]
Beta1 <- x[["Beta1"]]
Beta2 <- x[["Beta2"]]
Beta3 <- x[["Beta3"]]
Beta4 <- x[["Beta4"]]
BetaGreen <- x[["BetaGreen"]]
BetaMicro <- x[["BetaMicro"]]
BetaAna <- x[["BetaAna"]]
BetaNFix <- x[["BetaNFix"]]

sigma_p <- x[["sigma_p"]]

#Pull out environmental effects
Ntheta <- x[["Ntheta"]]
Ptheta <- x[["Ptheta"]]
Atheta <- x[["Atheta"]]
Dtheta <- x[["Dtheta"]]
Ttheta <- x[["Ttheta"]]
Ctheta <- x[["Ctheta"]]
Rtheta <- x[["Rtheta"]]

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 11

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  percentcover_latent[-c(1:2, 14:16, 29:33, 41:45), -1][1:time, ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  phos = stand_nut$oPhos_ug_P_L[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  amon = stand_nut$ammonium_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  discharge = discharge$stand_discharge[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  temp = stand_nut$temp_C[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  cond = stand_nut$cond_uS_cm[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  rad = swradiation$stand_rad[-c(1:2, 14:16, 29:33, 41:45)][1:time]
))

tox <- matrix(NA, runs, time)

# Build parameter matrixes
beta_matrix <- cbind(Beta0, Beta1, Beta2, Beta3, Beta4, Ntheta, Ptheta, Atheta, Dtheta, Ttheta,
                     Ctheta, Rtheta)
beta_lag_matrix <- cbind(BetaGreen, BetaMicro, BetaAna, BetaNFix)

for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  beta_lag <- beta_lag_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, beta%*%X1[t-1, ] + beta_lag%*%X1[t-2, 2:5],
                      sigma_p[z])
  }
}

pred2022 <- tox

sims2022median <- as.data.frame(apply(exp(pred2022)/1000, 2, median)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(time = 1:time) 
sims2022lquant <- as.data.frame(apply(exp(pred2022)/1000, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(time = 1:time)
sims2022uquant <- as.data.frame(apply(exp(pred2022)/1000, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(time = 1:time)

#Join together median and CI ranges
riversims2022 <- dplyr::left_join(sims2022median, sims2022lquant, by=c("time")) %>%
  dplyr::left_join(., sims2022uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 27, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7))


###
#Historical Predictions, 2023
###
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,12:13] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 12
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  percentcover_latent[-c(1:2, 14:16, 29:33, 41:45), -1][12:(11+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  amon = stand_nut$ammonium_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  discharge = discharge$stand_discharge[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  temp = stand_nut$temp_C[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  cond = stand_nut$cond_uS_cm[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  rad = swradiation$stand_rad[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)]
))


for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  beta_lag <- beta_lag_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, beta%*%X1[t-1, ] + beta_lag%*%X1[t-2, 2:5],
                      sigma_p[z])
  }
}

pred2023 <- tox

sims2023median <- as.data.frame(apply(exp(pred2023)/1000, 2, median)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(time = 1:time) 
sims2023lquant <- as.data.frame(apply(exp(pred2023)/1000, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(time = 1:time)
sims2023uquant <- as.data.frame(apply(exp(pred2023)/1000, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(time = 1:time)

#Join together median and CI ranges
riversims2023 <- dplyr::left_join(sims2023median, sims2023lquant, by=c("time")) %>%
  dplyr::left_join(., sims2023uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 27, year = 2023) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))



###
#Historical Predictions, 2024
###
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,24:30] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 7
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  percentcover_latent[, -1][24:(23+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  amon = stand_nut$ammonium_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  discharge = discharge$stand_discharge[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  temp = stand_nut$temp_C[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  cond = stand_nut$cond_uS_cm[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  rad = swradiation$stand_rad[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)]
))

for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  beta_lag <- beta_lag_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, beta%*%X1[t-1, ] + beta_lag%*%X1[t-2, 2:5],
                      sigma_p[z])
  }
}

pred2024 <- tox

sims2024median <- as.data.frame(apply(exp(pred2024)/1000, 2, median)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(time = 1:time) 
sims2024lquant <- as.data.frame(apply(exp(pred2024)/1000, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(time = 1:time)
sims2024uquant <- as.data.frame(apply(exp(pred2024)/1000, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(time = 1:time)

#Join together median and CI ranges
riversims2024 <- dplyr::left_join(sims2024median, sims2024lquant, by=c("time")) %>%
  dplyr::left_join(., sims2024uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 30, year = 2024) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))

#Join together simulation data
riverTACsimsallyears <- rbind(riversims2022, riversims2023, riversims2024)

###Create plot of River-wide predictions vs latent states

#Graphing palettes
#Create a color palette
mycols <- c("#791C55", "darkgreen")
mypal <- palette(mycols)
mypal <- palette(mycols)
names(mypal) = c("Latent", "Predicted")
rivercolScale <- scale_color_manual(name = "State Type", values = mypal)
filScale <- scale_fill_manual(name = "State Type", values = mypal)
linScale <- scale_linetype_manual(name = "State Type",
                                  values = c("Latent" = "11",
                                             "Predicted" = "solid"))

ggplot(riverTACsimsallyears, aes(x = model_date, y = toxins)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = "Predicted"), 
              alpha = 0.3) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", color = "Predicted"), size = 1.5) +
  # Latent points/lines
  geom_line(data = tox_params2_riverTAC,
            aes(y = median, linetype = "Latent", color = "Latent"), linewidth = 2) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(ylim = c(0,40)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "River-Wide Toxin") +
  rivercolScale + filScale + linScale + theme_bw()



#Compile model check dataframes into a single full timeseries matrix
predictives_TAC_toxins_river <- abind(pred2022, pred2023, pred2024, along = 2)
predictives_TAC_toxins_river <- exp(predictives_TAC_toxins_river)

#Save predictive output of Toxin River-Wide model
saveRDS(predictives_TAC_toxins_river, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Toxins_TAC_Pred_RiverWide.rds")) 


