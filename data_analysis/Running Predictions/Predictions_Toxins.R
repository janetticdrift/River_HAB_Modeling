#Historical Predictions of Anatoxin Concentrations, 2022
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(abind)

#Read in raw, uncleaned latent states and effect coefficients for toxin concentrations (River.fit and
  #Mat.fit), and latent states of percent cover and microscopy abundances 
  #(percentcover_latent and microscopyTM_latent)

                                    ###River-Wide###
#Percent cover predicted states
percentcover_models <- list(
  All = readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/River_predict_All.rds")),
  Biotic = readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/River_predict_Biotic.rds")),
  Abiotic = readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/River_predict_Abiotic.rds")),
  AbioticNoNut = readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/River_predict_Abioticnonut.rds")),
  TrueAbiotic = readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/River_predict_TrueAbiotic.rds"))
  )

#Microcoleus Anatoxin latent states
River.fit.TM_models <- list(
  All = readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_All_predictions.rds")),
  Biotic =readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_Biotic_predictions.rds")),
  Abiotic = readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_Abiotic_predictions.rds")),
  AbioticNoNut = readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_AbioticNoNut_predictions.rds")),
  TrueAbiotic = readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_TrueAbiotic_predictions.rds"))
  )

#Anabaena
River.fit.TAC <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TAC_River_predictions.rds"))

                                    ###Within-Mat###
microscopyTM_latent <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/WithinMat_Micro_LatentStates.rds"))
Mat.fit <- readRDS(here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_Mat_predictions.rds"))

#Read in cleaned toxin latent states, used for creating graphs (tox_params2)
source(here::here("data_analysis/Compare Obs Vs Modeled Outputs/Toxins_model_vs_real.R"))

######Simulate Toxins using microscopy data------

#Pull out community abundances and demographics 
x <- Mat.fit 
time <- 13
tox_conc <- x[["tox_raw"]][,1:2] #iterations, time
Beta0 <- x[["Beta0"]]
Beta1 <- x[["Beta1"]] #Anabaena
Beta2 <- x[["Beta2"]] #Epithemia
Beta3 <- x[["Beta3"]] #Geitlerinema

sigma_p <- x[["sigma_p"]]

Phi0 <- x[["Phi0"]]
PhiAna <- x[["PhiAna"]]

#Pull out environmental effects
Ntheta <- x[["Ntheta"]]
Ptheta <- x[["Ptheta"]]
Atheta <- x[["Atheta"]]
Dtheta <- x[["Dtheta"]]
Ttheta <- x[["Ttheta"]]
Ctheta <- x[["Ctheta"]]
Rtheta <- x[["Rtheta"]]

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  microscopyTM_latent[, -1][1:time, ],  #Abundances are on the log-scale
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)][1:time],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)][1:time],
  amon = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)][1:time],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)][1:time],
  temp = stand_nut$temp_C[-c(14:15, 29:30)][1:time],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)][1:time],
  rad = swradiation$stand_rad[-c(14:15, 29:30)][1:time]
))

#Inputs
runs <- length(Beta1) # number of model iterations
time <- 13

tox <- matrix(NA, runs, time)

# Build parameter matrixes
beta_matrix <- cbind(Beta0, Beta1, Beta2, Beta3, Ntheta, Ptheta, Atheta,
                     Dtheta, Ttheta,
                     Ctheta, Rtheta
)

for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    #----------
    #Hurdle: Do toxins initiate?
    #----------
    if(t <=2){ #Avoid indexing error with using Anabaena t-2 lag
      phi_t <- plogis(Phi0[z]) #plogis is the same inverse logit function as in Stan, just named differently
    } else {
      phi_t <- plogis(Phi0[z] + PhiAna[z] * X1[t-2,2])  #Use the Ana t-2 lag to predict initiation timing
    }
    
    #----------
    #Process: How many toxins are created?
    #----------
    #Probability of toxin initiation 
    toxin_initiate <- rbinom(1, size = 1, prob = phi_t) #Draw a random probability of toxin initiationg from the phi_t probability
    
    if(toxin_initiate == 1){ #If toxin production did initiate...
      tox[z,t] <- rnorm(1, mean = beta%*%X1[t-1, ], sd = sigma_p[z])
    } else {
      tox[z,t] <- log(0.0001) #If toxin production did not initiate...
    }
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
time <- 13
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  microscopyTM_latent[, -1][14:(13+time), ],  #Abundances are log-transformed
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
  
  #Simulation
  for (t in 3:time) {
    #----------
    #Hurdle: Do toxins initiate?
    #----------
    if(t <=2){ #Avoid indexing error with using Anabaena t-2 lag
      phi_t <- plogis(Phi0[z]) #plogis is the same inverse logit function as in Stan, just named differently
    } else {
      phi_t <- plogis(Phi0[z] + PhiAna[z] * X1[t-2,2])  #Use the Ana t-2 lag to predict initiation timing
    }
    
    #----------
    #Process: How many toxins are created?
    #----------
    #Probability of toxin initiation 
    toxin_initiate <- rbinom(1, size = 1, prob = phi_t) #Draw a random probability of toxin initiationg from the phi_t probability
    
    if(toxin_initiate == 1){ #If toxin production did initiate...
      tox[z,t] <- rnorm(1, mean = beta%*%X1[t-1, ], sd = sigma_p[z])
    } else {
      tox[z,t] <- log(0.0001) #If toxin production did not initiate...
    }
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
time <- 15
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  microscopyTM_latent[, -1][27:(26+time), ],  #Abundances are log-transformed
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
  
  #Simulation
  for (t in 3:time) {
    #----------
    #Hurdle: Do toxins initiate?
    #----------
    if(t <=2){ #Avoid indexing error with using Anabaena t-2 lag
      phi_t <- plogis(Phi0[z]) #plogis is the same inverse logit function as in Stan, just named differently
    } else {
      phi_t <- plogis(Phi0[z] + PhiAna[z] * X1[t-2,2])  #Use the Ana t-2 lag to predict initiation timing
    }
    
    #----------
    #Process: How many toxins are created?
    #----------
    #Probability of toxin initiation 
    toxin_initiate <- rbinom(1, size = 1, prob = phi_t) #Draw a random probability of toxin initiationg from the phi_t probability
    
    if(toxin_initiate == 1){ #If toxin production did initiate...
      tox[z,t] <- rnorm(1, mean = beta%*%X1[t-1, ], sd = sigma_p[z])
    } else {
      tox[z,t] <- log(0.0001) #If toxin production did not initiate...
    }
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
matsimsallyears <- rbind(matsims2022, matsims2023, matsims2024) %>% 
  left_join(tox_params2_mat %>% dplyr::select(model_date, median),
            by = "model_date") %>% 
  dplyr::rename(Predicted = toxins, Latent = median) %>% 
  pivot_longer(cols = c(Predicted, Latent), names_to = "StateType",
               values_to = "toxins")

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

ggplot(matsimsallyears, aes(x = model_date, y = toxins, color = StateType,
                                fill = StateType, linetype = StateType)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(data = filter(matsimsallyears, StateType == "Predicted"), 
              aes(ymin = `CIlower`, ymax = `CIupper`), color = NA,
              alpha = 0.3) +
  # Predicted points/lines
  geom_line(data = filter(matsimsallyears, StateType == "Predicted"), size = 1.5) +
  # Latent points/lines
  geom_line(data = filter(matsimsallyears, StateType == "Latent"), linewidth = 2) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(ylim = c(0,80)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "Within-Mat: Latent vs. Predicted Toxin Concentrations from Microcoleus Mats") +
  matcolScale + filScale + linScale + theme_bw()



#Compile model check dataframes into a single full timeseries matrix
predictives_toxins_mats1 <- abind(pred2022, pred2023, pred2024, along = 2)
predictives_toxins_mats <- exp(predictives_toxins_mats1)

#Save predictive output of Toxin Within-Mat model
saveRDS(predictives_toxins_mats, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Toxins_Pred_Withinmat.rds")) 




######Simulate Microcoleus Toxins using percent cover data-----

#Create vector of model names
model_names <- c("All","Biotic","Abiotic","AbioticNoNut", "TrueAbiotic")
#Create empty list to store the completed simulation results for each model
all_model_results <- list()
#Create empty list to store the plots for each model
all_model_plots <- list()

#Save the timesteps that should be removed from environmnetal and algal abundance data
env_keep <- -c(14:15, 29:30)

#Index out each environmental variable
nitrate_clean <- stand_nut$nitrate_mg_N_L[env_keep]
phos_clean <- stand_nut$oPhos_ug_P_L[env_keep]
amon_clean <- stand_nut$ammonium_mg_N_L[env_keep]
discharge_clean <- discharge$stand_discharge[env_keep]
temp_clean <- stand_nut$temp_C[env_keep]
cond_clean <- stand_nut$cond_uS_cm[env_keep]
rad_clean <- swradiation$stand_rad[env_keep]

################################################
#Use the simulate_toxin_year function to run simulations for All, Biotic, 
#Abiotic, Abiotic No Nutrients, and True Abiotic models
################################################

#Read in the simulate_toxin_year function
source(here::here("/data cleaning/Functions.R"))

#Use a for loop to run through the 3 years of simulations for each river-wide model
for (model in model_names) {
  #Extract the model
  x <- River.fit.TM_models[[model]]
  #Extract the corresponding predicted state of percent cover abundances
  percentcover_predicted <- percentcover_models[[model]]
  #Extract the corresponding latent state of toxins, from Toxins_model_vs_real.R, for plotting
  tox_params2_index <- tox_params2[[model]]
  
  #2022 simulations
  result2022 <- simulate_toxin_year(x = x,
                                    percentcover_predicted = percentcover_predicted,
                                    year = 2022,
                                    start_timestep = 1,
                                    time = 13,
                                    week_offset = 25)

  #2023 simulations
  result2023 <- simulate_toxin_year(x = x, 
                                    percentcover_predicted = percentcover_predicted,
                                    year = 2023,
                                    start_timestep = 14,
                                    time = 13,
                                    week_offset = 26)

  #2024 simulations
  result2024 <- simulate_toxin_year(x = x,
                                    percentcover_predicted = percentcover_predicted,
                                    year = 2024,
                                    start_timestep = 27,
                                    time = 15,
                                    week_offset = 26)
  
  #Pull out the summarized predictions for each year that function was used on
  riversims2022 <- result2022$predictions
  riversims2023 <- result2023$predictions
  riversims2024 <- result2024$predictions
  
  #Pull out the raw predictions for each year that the function was used on; used
    #for calculating model fit indices
  pred2022 <- result2022$raw
  pred2023 <- result2023$raw
  pred2024 <- result2024$raw
  
  #Merge summarized predictions from all years into one dataframe, and then also merge
    #with dataframe containing true latent state values
  riverTMsimsallyears <- rbind(riversims2022, riversims2023, riversims2024) %>% 
    dplyr::left_join(tox_params2_index %>% 
                       dplyr::select(model_date, median), by = "model_date") %>% 
    dplyr::rename(Predicted = toxins, Latent = median) %>% 
    tidyr::pivot_longer(cols = c(Predicted, Latent), names_to = "StateType",values_to = "toxins")
  
  #Store raw results in an empty list
  all_model_results[[model]] <- list(simulations = riverTMsimsallyears,
                                     pred2022 = pred2022,
                                     pred2023 = pred2023,
                                     pred2024 = pred2024)
 
   
                                    ######Create Figures#####
  
  #Color palette
  mycols <- c("#791C55","#41789A")
  mypal <- palette(mycols)
  mypal <- palette(mycols)
  names(mypal) <- c("Latent","Predicted")
  riverTMcolScale <- scale_color_manual(name = "State Type", values = mypal)
  filScale <- scale_fill_manual(name = "State Type", values = mypal)
  linScale <- scale_linetype_manual(name = "State Type", values = c("Latent" = "11",
                                                                    "Predicted" = "solid"))
  
  
  p <- ggplot(riverTMsimsallyears, aes(x = model_date, y = toxins, color = StateType, 
                                       fill = StateType, linetype = StateType)) +
    facet_wrap(~year, scales = "free_x") +
    #Predicted CI
    geom_ribbon(data = filter(riverTMsimsallyears, StateType == "Predicted"),
                aes(ymin = CIlower, ymax = CIupper), color = NA, alpha = 0.3) +
    #Predicted median
    geom_line(data = filter(riverTMsimsallyears, StateType == "Predicted"),
              linewidth = 1.5) +
    #Latent state line
    geom_line(data = filter(riverTMsimsallyears, StateType == "Latent"), linewidth = 2) +
    scale_y_continuous(breaks = seq(0, 100, 10)) +
    coord_cartesian(ylim = c(0, 80)) +
    #Unique title per plot
    labs(x = "Date", y = "Anatoxin Concentration (ug/g)", 
         title = paste0("River-Wide TM Mats: Latent vs. Predicted Toxin Concentrations - ",
                        model," Model")) +
    riverTMcolScale + filScale + linScale + theme_bw()
  
  #Display figure
  print(p)
  
  #Store figure in list
  all_model_plots[[model]] <- p
  
  #Compile raw predictions for model fit checks
  
  predictives_TM_toxins_river1 <- abind(pred2022, pred2023, pred2024, along = 2)
  predictives_TM_toxins_river <- exp(predictives_TM_toxins_river1) #Back-transform from log scale
  
  #Save raw prediction files for model fit checks
  predictive_file <- here::here(paste0("data/Outputs for Sims and Model Fits/Predicted States/",
                                       "Toxins_TM_Pred_RiverWide_", model,".rds"))
  
  saveRDS(predictives_TM_toxins_river, file = predictive_file)

}


######Simulate Anabaena Toxins using percent cover data-----

#Pull out community abundances and demographics 
x <- River.fit.TAC 
tox_conc <- x[["tox_raw"]][,1:2] #iterations, time
Beta0 <- x[["Beta0"]]
Beta1 <- x[["Beta1"]]
Beta2 <- x[["Beta2"]]
Beta3 <- x[["Beta3"]]
Beta4 <- x[["Beta4"]]

sigma_p <- x[["sigma_p"]]

Phi0 <- x[["Phi0"]]
PhiAna <- x[["PhiAna"]]

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
percentcover_predicted <- percentcover_models[["All"]]

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  percentcover_predicted[-c(1:2, 14:16, 29:33, 41:45), -1][1:time, ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  phos = stand_nut$oPhos_ug_P_L[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  amon = stand_nut$ammonium_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  discharge = discharge$stand_discharge[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  temp = stand_nut$temp_C[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  cond = stand_nut$cond_uS_cm[-c(1:2, 14:16, 29:33, 41:45)][1:time],
  rad = swradiation$stand_rad[-c(1:2, 14:16, 29:33, 41:45)][1:time]
  # DIN = stand_nut$DIN[-c(1:2, 14:16, 29:33, 41:45)][1:time]
))

tox <- matrix(NA, runs, time)

# Build parameter matrixes
beta_matrix <- cbind(Beta0, Beta1, Beta2, Beta3, Beta4, Ntheta, Ptheta, Atheta,
                     Dtheta, Ttheta,
                     Ctheta, Rtheta
                     )

for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    #----------
    #Hurdle: Do toxins initiate?
    #----------
    if(t <=2){ #Avoid indexing error with using Anabaena t-2 lag
      phi_t <- plogis(Phi0[z]) #plogis is the same inverse logit function as in Stan, just named differently
    } else {
      phi_t <- plogis(Phi0[z] + PhiAna[z] * X1[t-2,4])  #Use the Ana t-2 lag to predict initiation timing
    }
    
    #----------
    #Process: How many toxins are created?
    #----------
    #Probability of toxin initiation 
    toxin_initiate <- rbinom(1, size = 1, prob = phi_t) #Draw a random probability of toxin initiationg from the phi_t probability
    
    if(toxin_initiate == 1){ #If toxin production did initiate...
      tox[z,t] <- rnorm(1, mean = beta%*%X1[t-1, ], sd = sigma_p[z])
    } else {
      tox[z,t] <- log(0.0001) #If toxin production did not initiate...
    }
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
time <- 12
tox <- matrix(NA, runs, time)

#Create design matrix
X1 <- as.matrix(cbind(
  intercept = rep(1, time),
  percentcover_predicted[-c(1:2, 14:16, 29:33, 41:45), -1][12:(11+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  amon = stand_nut$ammonium_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  discharge = discharge$stand_discharge[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  temp = stand_nut$temp_C[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  cond = stand_nut$cond_uS_cm[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)],
  rad = swradiation$stand_rad[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)]
  # DIN = stand_nut$DIN[-c(1:2, 14:16, 29:33, 41:45)][12:(11+time)]
))


for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    #----------
    #Hurdle: Do toxins initiate?
    #----------
    if(t <=2){ #Avoid indexing error with using Anabaena t-2 lag
      phi_t <- plogis(Phi0[z]) #plogis is the same inverse logit function as in Stan, just named differently
    } else {
      phi_t <- plogis(Phi0[z] + PhiAna[z] * X1[t-2,4])  #Use the Ana t-2 lag to predict initiation timing
    }
    
    #----------
    #Process: How many toxins are created?
    #----------
    #Probability of toxin initiation 
    toxin_initiate <- rbinom(1, size = 1, prob = phi_t) #Draw a random probability of toxin initiationg from the phi_t probability
    
    if(toxin_initiate == 1){ #If toxin production did initiate...
      tox[z,t] <- rnorm(1, mean = beta%*%X1[t-1, ], sd = sigma_p[z])
    } else {
      tox[z,t] <- log(0.0001) #If toxin production did not initiate...
    }
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
  percentcover_predicted[, -1][24:(23+time), ],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  phos = stand_nut$oPhos_ug_P_L[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  amon = stand_nut$ammonium_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  discharge = discharge$stand_discharge[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  temp = stand_nut$temp_C[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  cond = stand_nut$cond_uS_cm[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)],
  rad = swradiation$stand_rad[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)]
  # DIN = stand_nut$DIN[-c(1:2, 14:16, 29:33, 41:45)][24:(23+time)]
))

for (z in 1:runs) {
  
  #Set initial tox concentrations for the first two skipped days
  tox[z,1] <- log(tox_conc[z,1] + 1e-6)
  tox[z,2] <- log(tox_conc[z,2] + 1e-6)
  
  #Identify timesteps in coefficient matrices
  beta <- beta_matrix[z, ]
  
  #Simulation
  for (t in 3:time) {
    #----------
    #Hurdle: Do toxins initiate?
    #----------
    if(t <=2){ #Avoid indexing error with using Anabaena t-2 lag
      phi_t <- plogis(Phi0[z]) #plogis is the same inverse logit function as in Stan, just named differently
    } else {
      phi_t <- plogis(Phi0[z] + PhiAna[z] * X1[t-2,4])  #Use the Ana t-2 lag to predict initiation timing
    }
    
    #----------
    #Process: How many toxins are created?
    #----------
    #Probability of toxin initiation 
    toxin_initiate <- rbinom(1, size = 1, prob = phi_t) #Draw a random probability of toxin initiationg from the phi_t probability
    
    if(toxin_initiate == 1){ #If toxin production did initiate...
      tox[z,t] <- rnorm(1, mean = beta%*%X1[t-1, ], sd = sigma_p[z])
    } else {
      tox[z,t] <- log(0.0001) #If toxin production did not initiate...
    }
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
riverTACsimsallyears <- rbind(riversims2022, riversims2023, riversims2024) %>% 
  left_join(tox_params2_riverTAC %>% dplyr::select(model_date, median),
            by = "model_date") %>% 
  dplyr::rename(Predicted = toxins, Latent = median) %>% 
  pivot_longer(cols = c(Predicted, Latent), names_to = "StateType",
               values_to = "toxins")

###Create plot of River-wide predictions vs latent states

#Graphing palettes
#Create a color palette
mycols <- c("#791C55", "darkgreen")
mypal <- palette(mycols)
mypal <- palette(mycols)
names(mypal) = c("Latent", "Predicted")
riverTACcolScale <- scale_color_manual(name = "State Type", values = mypal)
filScale <- scale_fill_manual(name = "State Type", values = mypal)
linScale <- scale_linetype_manual(name = "State Type",
                                  values = c("Latent" = "11",
                                             "Predicted" = "solid"))

ggplot(riverTACsimsallyears, aes(x = model_date, y = toxins, color = StateType,
                                 fill = StateType, linetype = StateType)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(data = filter(riverTACsimsallyears, StateType == "Predicted"), 
              aes(ymin = `CIlower`, ymax = `CIupper`), color = NA,
              alpha = 0.3) +
  # Predicted points/lines
  geom_line(data = filter(riverTACsimsallyears, StateType == "Predicted"), size = 1.5) +
  # Latent points/lines
  geom_line(data = filter(riverTACsimsallyears, StateType == "Latent"), linewidth = 2) +
  scale_y_continuous(breaks = seq(0, 150, 10)) +
  coord_cartesian(ylim = c(0,150)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "River-Wide: Latent vs. Predicted Toxin Concentrations from Anabaena Mats") +
  riverTACcolScale + filScale + linScale + theme_bw()



#Compile model check dataframes into a single full timeseries matrix
predictives_TAC_toxins_river <- abind(pred2022, pred2023, pred2024, along = 2)
predictives_TAC_toxins_river <- exp(predictives_TAC_toxins_river)

#Save predictive output of Toxin River-Wide model
saveRDS(predictives_TAC_toxins_river, 
        file = here::here("data/Outputs for Sims and Model Fits/Predicted States/Toxins_TAC_Pred_RiverWide.rds")) 


