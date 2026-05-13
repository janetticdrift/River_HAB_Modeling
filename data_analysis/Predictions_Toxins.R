#Historical Predictions of Anatoxin Concentrations, 2022
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(abind)

#Read in latent states and effect coefficients
River.fit <- readRDS(here::here("data/Model Fits/Anatoxin_River_predictions.rds"))
Mat.fit <- readRDS(here::here("data/Model Fits/Anatoxin_Mat_predictions.rds"))


######Simulate Toxins using microscopy data

#Pull out community abundances and demographics 
x <- Mat.fit 
tox_conc <- x[["tox_raw"]][,1:2] #iterations, time
Beta0 <- x[["Beta0"]]
Beta1 <- x[["Beta1"]]
Beta2 <- x[["Beta2"]]
Beta3 <- x[["Beta3"]]
BetaAna  <- x[["BetaAna"]]
BetaEpi  <- x[["BetaEpi"]]
BetaGeit <- x[["BetaGeit"]]

sigma_p <- x[["sigma_p"]]

#Pull out environmental effects
Ntheta <- x[["Ntheta"]]

Ptheta <- x[["Ptheta"]]

Atheta <- x[["Atheta"]]

Dtheta <- x[["Dtheta"]]

Ttheta <- x[["Ttheta"]]

Ctheta <- x[["Ctheta"]]

Rtheta <- x[["Rtheta"]]

#Create design matrix
X1 <- cbind(
  intercept = rep(1, time),
  TM_latent[, -1][1:time, ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][1:time],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][1:time],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][1:time],
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
    Dtheta[z],
    Ttheta[z],
    Ctheta[z],
    Rtheta[z]
  )
  
  beta_lag <- c(
    BetaAna[z],
    BetaEpi[z],
    BetaGeit[z]
  )
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, sum(X1[t-1, ] * beta) +
                        sum(X1[t-2, 2:4] * beta_lag)
                      , sigma_p[z])
  }
}

pred2022 <- tox

sims2022mean <- as.data.frame(apply(tox, 2, mean)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(toxins = exp(toxins)) %>% 
  dplyr::mutate(time = 1:time) 
sims2022lquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(CIlower = exp(CIlower)) %>% 
  dplyr::mutate(time = 1:time)
sims2022uquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(CIupper = exp(CIupper)) %>% 
  dplyr::mutate(time = 1:time)

#Join together mean and CI ranges
matsims2022 <- dplyr::left_join(sims2022mean, sims2022lquant, by=c("time")) %>%
  dplyr::left_join(., sims2022uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))


#####################################
#Historical Predictions, 2023
#####################################
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,14:15] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 13
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- cbind(
  intercept = rep(1, time),
  TM_latent[, -1][14:(13+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][14:(13+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][14:(13+time)],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][14:(13+time)],
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
    Dtheta[z],
    Ttheta[z],
    Ctheta[z],
    Rtheta[z]
  )
  
  beta_lag <- c(
    BetaAna[z],
    BetaEpi[z],
    BetaGeit[z]
  )
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, sum(X1[t-1, ] * beta) +
                        sum(X1[t-2, 2:4] * beta_lag)
                      , sigma_p[z])
  }
}

pred2023 <- tox

sims2023mean <- as.data.frame(apply(tox, 2, mean)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(toxins = exp(toxins)) %>% 
  dplyr::mutate(time = 1:time) 
sims2023lquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(CIlower = exp(CIlower)) %>% 
  dplyr::mutate(time = 1:time)
sims2023uquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(CIupper = exp(CIupper)) %>% 
  dplyr::mutate(time = 1:time)

#Join together mean and CI ranges
matsims2023 <- dplyr::left_join(sims2023mean, sims2023lquant, by=c("time")) %>%
  dplyr::left_join(., sims2023uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 26, year = 2023) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))



#####################################
#Historical Predictions, 2024
#####################################
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,27:28] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 15
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- cbind(
  intercept = rep(1, time),
  TM_latent[, -1][27:(26+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][27:(26+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][27:(26+time)],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][27:(26+time)],
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
    Dtheta[z],
    Ttheta[z],
    Ctheta[z],
    Rtheta[z]
  )
  
  beta_lag <- c(
    BetaAna[z],
    BetaEpi[z],
    BetaGeit[z]
  )
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, sum(X1[t-1, ] * beta) +
                        sum(X1[t-2, 2:4] * beta_lag)
                      , sigma_p[z])
  }
}

pred2024 <- tox

sims2024mean <- as.data.frame(apply(tox, 2, mean)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(toxins = exp(toxins)) %>% 
  dplyr::mutate(time = 1:time) 
sims2024lquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(CIlower = exp(CIlower)) %>% 
  dplyr::mutate(time = 1:time)
sims2024uquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(CIupper = exp(CIupper)) %>% 
  dplyr::mutate(time = 1:time)

#Join together mean and CI ranges
matsims2024 <- dplyr::left_join(sims2024mean, sims2024lquant, by=c("time")) %>%
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
mycols <- c("lightsalmon2", "brown")
mypal <- palette(mycols)
names(mypal) = c("Latent", "Predicted")
colScale <- scale_color_manual(name = "State Type", values = mypal)
filScale <- scale_fill_manual(name = "State Type", values = mypal)
linScale <- scale_linetype_manual(name = "State Type",
                                  values = c("Latent" = "11",
                                             "Predicted" = "solid"))

ggplot(matsimsallyears, aes(x = model_date, y = toxins)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = "Predicted"), 
              alpha = 0.3) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", color = "Predicted"), size = 1.5) +
  # Latent points/lines
  geom_line(data = tox_params2_mat,
            aes(y = mean, linetype = "Latent", color = "Latent"), linewidth = 2) +
  scale_y_continuous(breaks = seq(0, 200, 10)) +
  coord_cartesian(ylim = c(0,160)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "Within-Mat: Latent vs. Predicted Toxin Concentrations") +
  colScale + filScale + linScale



#Compile model check dataframes into a single full timeseries matrix
predictives_toxins_mats <- abind(pred2022, pred2023, pred2024, along = 2)
predictives_toxins_mats <- exp(predictives_toxins_mats)

#Save predictive output of Toxin Within-Mat model
saveRDS(predictives_toxins_mats, 
        file = here::here("data/Toxins_Pred_Withinmat.rds"))




######Simulate Toxins using percent cover data

#Pull out community abundances and demographics 
x <- River.fit 
tox_conc <- x[["tox_raw"]][,1:2] #iterations, time
Beta0 <- x[["Beta0"]]
Beta1 <- x[["Beta1"]]
Beta2 <- x[["Beta2"]]
Beta3 <- x[["Beta3"]]
Beta4 <- x[["Beta4"]]
BetaGreen  <- x[["BetaGreen"]]
BetaMicro  <- x[["BetaMicro"]]
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
X1 <- cbind(
  intercept = rep(1, time),
  River_latent[-c(14:15, 29:30), -1][1:time, ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][1:time],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][1:time],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][1:time],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)][1:time],
  temp = stand_nut$temp_C[-c(14:15, 29:30)][1:time],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)][1:time],
  rad = swradiation$stand_rad[-c(14:15, 29:30)][1:time]
)

tox <- matrix(NA, runs, time)

for (z in 1:runs) {
  
  # Build parameter vectors
  beta <- c(
    Beta0[z],
    Beta1[z],
    Beta2[z],
    Beta3[z],
    Beta4[z],
    Ntheta[z],
    Ptheta[z],
    Atheta[z],
    Dtheta[z],
    Ttheta[z],
    Ctheta[z],
    Rtheta[z]
  )
  
  beta_lag <- c(
    BetaGreen[z],
    BetaMicro[z],
    BetaAna[z],
    BetaNFix[z]
  )
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, sum(X1[t-1, ] * beta) +
                        sum(X1[t-2, 2:4] * beta_lag)
                      , sigma_p[z])
  }
}

pred2022 <- tox

sims2022mean <- as.data.frame(apply(tox, 2, mean)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(toxins = exp(toxins)) %>% 
  dplyr::mutate(time = 1:time) 
sims2022lquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(CIlower = exp(CIlower)) %>% 
  dplyr::mutate(time = 1:time)
sims2022uquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(CIupper = exp(CIupper)) %>% 
  dplyr::mutate(time = 1:time)

#Join together mean and CI ranges
matsims2022 <- dplyr::left_join(sims2022mean, sims2022lquant, by=c("time")) %>%
  dplyr::left_join(., sims2022uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 25, year = 2022) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))


#####################################
#Historical Predictions, 2023
#####################################
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,14:15] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 13
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- cbind(
  intercept = rep(1, time),
  TM_latent[, -1][14:(13+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][14:(13+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][14:(13+time)],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][14:(13+time)],
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
    Dtheta[z],
    Ttheta[z],
    Ctheta[z],
    Rtheta[z]
  )
  
  beta_lag <- c(
    BetaAna[z],
    BetaEpi[z],
    BetaGeit[z]
  )
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, sum(X1[t-1, ] * beta) +
                        sum(X1[t-2, 2:4] * beta_lag)
                      , sigma_p[z])
  }
}

pred2023 <- tox

sims2023mean <- as.data.frame(apply(tox, 2, mean)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(toxins = exp(toxins)) %>% 
  dplyr::mutate(time = 1:time) 
sims2023lquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(CIlower = exp(CIlower)) %>% 
  dplyr::mutate(time = 1:time)
sims2023uquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(CIupper = exp(CIupper)) %>% 
  dplyr::mutate(time = 1:time)

#Join together mean and CI ranges
matsims2023 <- dplyr::left_join(sims2023mean, sims2023lquant, by=c("time")) %>%
  dplyr::left_join(., sims2023uquant, by=c("time")) %>% 
  dplyr::mutate(real_week = time + 26, year = 2023) %>% 
  dplyr::mutate(model_date = ceiling_date(
    ymd(paste(year, "01", "01", sep = "-")) + (real_week - 1) * 7 - 1,
    "week", week_start = 7
  ))



#####################################
#Historical Predictions, 2024
#####################################
#Pull out new starting year
tox_conc <- x[["tox_raw"]][,27:28] #iterations, time

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 15
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- cbind(
  intercept = rep(1, time),
  TM_latent[, -1][27:(26+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][27:(26+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][27:(26+time)],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][27:(26+time)],
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
    Dtheta[z],
    Ttheta[z],
    Ctheta[z],
    Rtheta[z]
  )
  
  beta_lag <- c(
    BetaAna[z],
    BetaEpi[z],
    BetaGeit[z]
  )
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Simulation
  for (t in 3:time) {
    
    tox[z,t] <- rnorm(1, sum(X1[t-1, ] * beta) +
                        sum(X1[t-2, 2:4] * beta_lag)
                      , sigma_p[z])
  }
}

pred2024 <- tox

sims2024mean <- as.data.frame(apply(tox, 2, mean)) %>% 
  dplyr::rename(toxins = 1) %>% 
  dplyr::mutate(toxins = exp(toxins)) %>% 
  dplyr::mutate(time = 1:time) 
sims2024lquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.025)) %>% 
  dplyr::rename(CIlower = 1) %>% 
  dplyr::mutate(CIlower = exp(CIlower)) %>% 
  dplyr::mutate(time = 1:time)
sims2024uquant <- as.data.frame(apply(tox, 2, quantile, probs = 0.975)) %>% 
  dplyr::rename(CIupper = 1) %>% 
  dplyr::mutate(CIupper = exp(CIupper)) %>% 
  dplyr::mutate(time = 1:time)

#Join together mean and CI ranges
matsims2024 <- dplyr::left_join(sims2024mean, sims2024lquant, by=c("time")) %>%
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
mycols <- c("lightsalmon2", "brown")
mypal <- palette(mycols)
names(mypal) = c("Latent", "Predicted")
colScale <- scale_color_manual(name = "State Type", values = mypal)
filScale <- scale_fill_manual(name = "State Type", values = mypal)
linScale <- scale_linetype_manual(name = "State Type",
                                  values = c("Latent" = "11",
                                             "Predicted" = "solid"))

ggplot(matsimsallyears, aes(x = model_date, y = toxins)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = "Predicted"), 
              alpha = 0.3) +
  # Predicted points/lines
  geom_line(aes(linetype = "Predicted", color = "Predicted"), size = 1.5) +
  # Latent points/lines
  geom_line(data = tox_params2_mat,
            aes(y = mean, linetype = "Latent", color = "Latent"), linewidth = 2) +
  scale_y_continuous(breaks = seq(0, 200, 10)) +
  coord_cartesian(ylim = c(0,160)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "Within-Mat: Latent vs. Predicted Toxin Concentrations") +
  colScale + filScale + linScale



#Compile model check dataframes into a single full timeseries matrix
predictives_toxins_mats <- abind(pred2022, pred2023, pred2024, along = 2)
predictives_toxins_mats <- exp(predictives_toxins_mats)

#Save predictive output of Toxin Within-Mat model
saveRDS(predictives_toxins_mats, 
        file = here::here("data/Toxins_Pred_Withinmat.rds"))
