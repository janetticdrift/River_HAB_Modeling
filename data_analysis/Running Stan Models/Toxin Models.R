########################################
###Anatoxin Dataframes for STAN Model###
########################################
#This file contains code reading in and cleaning data on anatoxins, gene count, microscopy data,
#and environmental covariates to predict anatoxins

#Packages----
library(rstan)
library(tidyverse)
library(dataRetrieval)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(loo)

#Read in environmental, microscopy, and anatoxin data
source(here::here("data_cleaning/cleaning_HAB.R"))

#Clean dataframes to feed into the toxin model

#Isolate Microcoleus mats
toxins <- toxindf %>% 
  dplyr::filter(sample_type == "Microcoleus") %>% 
  dplyr::mutate(field_date = replace(field_date, field_date == as.Date("2023-07-11"), 
                                     as.Date("2023-07-10"))) %>%  #Replace 2023-07-11 with 07/10 so they are on the same week
  dplyr::mutate(field_date = replace(field_date, field_date == as.Date("2022-09-06"),
                                     as.Date("2022-09-08"))) %>%  #Replace 2022/09/06 to 09/08 so they are on the same week 
  dplyr::filter(field_date != as.Date("2024-06-19")) %>% #Remove 6/19/24 as it wasn't sampled in microscopy data
  dplyr::mutate(field_date = as.Date(field_date, format="%m/%d/%y"))

#Save cleaned output for visualizing observational vs latent states
saveRDS(toxins, 
        file = here::here("data/Outputs for Obs vs Real/obs_toxins.rds"))


anaCsplit <- toxins %>% 
  ungroup() %>% 
  group_split(year)
year1 <- anaCsplit[[1]]
year2 <- anaCsplit[[2]]
year3 <- anaCsplit[[3]]

#year 2022
year1_index <- year1 %>% 
  dplyr::mutate(timestep = dense_rank(field_date)) %>% 
  dplyr::mutate(week = case_when(timestep == 1 ~ 1,
                                 timestep == 2 ~ 3,
                                 timestep == 3 ~ 5,
                                 timestep == 4 ~ 7,
                                 timestep == 5 ~ 9,
                                 timestep == 6 ~ 11,
                                 timestep == 7 ~ 13))

#year 2023 - Note that they did not sample anatoxins the first week of %covers
year2_index <- year2 %>% 
  dplyr::mutate(timestep = dense_rank(field_date)) %>% 
  dplyr::mutate(week = timestep)

#year 2024
year3_index <- year3 %>% 
  dplyr::mutate(timestep = dense_rank(field_date)) %>%
  dplyr::mutate(week = case_when(timestep == 1 ~ 1,
                                 timestep == 2 ~ 3,
                                 timestep == 3 ~ 5,
                                 timestep == 4 ~ 7,
                                 timestep == 5 ~ 9,
                                 timestep == 6 ~ 11,
                                 timestep == 7 ~ 13,
                                 timestep == 8 ~ 15,
                                 timestep == 9 ~ 17))

atx <- rbind(year1_index, year2_index, year3_index) %>% 
  dplyr::select(-c(field_date, sample_type, timestep)) %>% 
  group_by(year) %>% 
  complete(nesting(reach), week = seq(min(week), max(week), 1L)) %>% 
  ungroup() %>%
  dplyr::mutate(reach = as.numeric(factor(reach))) 

#Gather data into stan list format
anatoxin_data <- atx %>% 
  group_by(year, week) %>% 
  dplyr::summarise(ATX_all_ug_g = mean(ATX_all_ug_g, na.rm = TRUE)) %>% #Average across reaches, removing reaches where no ATX was collected
  dplyr::mutate(ATX_all_ug_g = round(ATX_all_ug_g, digits = 3), #Editing data for poisson distribution
                ATX_all_ug_g = ATX_all_ug_g*1000) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) %>% 
  mutate(firstday = if_else(week == 1 & (year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) %>% 
  dplyr::mutate(is_obs  = ifelse(ATX_all_ug_g == -99, 0, 1), #Editing data for poisson
                ATX_all_ug_g = ifelse(ATX_all_ug_g == -99, 0, ATX_all_ug_g))

#Save cleaned output for visualizing observational vs latent states
saveRDS(anatoxin_data, 
        file = here::here("data/Outputs for Obs vs Real/stanformat_toxins.rds"))

#---------------------------------------------------------------------------------------
#CREATE MODEL FOR MICROCOLEUS WITHIN-MAT DATA

#Gather latent states of microscopy abundances from Within-Mat model
matmodel_TM <- readRDS(here::here("data/Outputs for Obs vs Real/WithinMat_Micro.rds"))

TM_latent1 <- as.data.frame(matmodel_TM) %>% 
  dplyr::select(matches("n\\[")) %>% 
  t 

TM_latent <- as.data.frame(TM_latent1) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE)) %>% 
  dplyr::mutate(Species = case_when(grepl("[1,", group, fixed=TRUE) ~ 'Anabaena',
                                    grepl("[2,", group, fixed=TRUE) ~ 'Epithemia Diatoms',
                                    grepl("[3,", group, fixed=TRUE) ~ 'Geitlerinema')) %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,2])) %>% 
  dplyr::select(-group) %>% 
  pivot_wider(names_from = Species, values_from = median) %>% 
  arrange(time)

#Save cleaned microscopy latent dataframe for making toxin prediction simulation
saveRDS(TM_latent, 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/WithinMat_Micro_LatentStates.rds"))

#Create design matrix
X1 <- cbind(
  intercept = 1,
  TM_latent[, -1],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)],
  temp = stand_nut$temp_C[-c(14:15, 29:30)],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)],
  rad = swradiation$stand_rad[-c(14:15, 29:30)]
)

#Combine with other information into model list
model.atx.mat <- list("uniqueID" = nrow(anatoxin_data),
                  "is_obs" = anatoxin_data$is_obs, #poisson edit
                  "firstdays" = anatoxin_data$firstday,
                  "Toxins" = as.integer(anatoxin_data$ATX_all_ug_g), #poisson edit needs as.integer
                  "Nspecies" = as.integer(ncol(TM_latent)-1),
                  "X1" = X1,
                  "Npredictors" = ncol(X1)
)

#---------------------------------------------------------------------------------------
#River-Wide 
#Gather latent states of percent cover abundances
rivermodel <- readRDS(here::here("data/Outputs for Obs vs Real/Riverwide_AllVariables.rds"))

River_latent1 <- as.data.frame(rivermodel) %>% 
  dplyr::select(matches("n\\[")) %>% 
  t 

River_latent <- as.data.frame(River_latent1) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE)) %>% 
  dplyr::mutate(Species = case_when(grepl("[1,", group, fixed=TRUE) ~ 'Green Algae',
                                    grepl("[2,", group, fixed=TRUE) ~ 'Microcoleus',
                                    grepl("[3,", group, fixed=TRUE) ~ 'Anabaena',
                                    grepl("[4,", group, fixed=TRUE) ~ 'Other N Fixers')) %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,2])) %>% 
  dplyr::select(-group) %>% 
  pivot_wider(names_from = Species, values_from = median) %>% 
  arrange(time)

#Save cleaned percent cover latent dataframe for making toxin prediction simulation
saveRDS(River_latent, 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_LatentStates.rds"))

#Create design matrix
X2 <- cbind(
  intercept = 1,
  River_latent[-c(14:15, 29:30), -1],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)],
  temp = stand_nut$temp_C[-c(14:15, 29:30)],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)],
  rad = swradiation$stand_rad[-c(14:15, 29:30)]
)

#Combine with other information into model list
model.atx.river <- list("uniqueID" = nrow(anatoxin_data),
                  "is_obs" = anatoxin_data$is_obs, #poisson edit, was this an observed day?
                  "firstdays" = anatoxin_data$firstday,
                  "Toxins" = as.integer(anatoxin_data$ATX_all_ug_g), #poisson edit needs as.integer
                  "Nspecies" = as.integer(ncol(River_latent)-1),
                  "X2" = X2,
                  "Npredictors" = ncol(X2)
)

#-------------------------------------------------------------------------------------------------
#Run models

setwd(here::here("data_analysis/Stan Models")) #Set working directory to current folder

options(mc.cores = parallel::detectCores())

#Set starting values
init_fun_atx <- function() list(
  sigma_p = 0.5,    
  Beta0 = 0,
  Beta1 = 0,     # small start for species abundances
  Beta2 = 0,
  Beta3 = 0,
  tox_nc = rep(0, nrow(anatoxin_data)) #nrow is the time length
 )

#Estimate anatoxins in river-wide assemblages
fit.atx.river <-  stan(file = "HAB_toxins_River_Wide.stan", data = model.atx.river, chains = 3, iter = 6000,
                       warmup = 3000, refresh=100, init = init_fun_atx, control = list(adapt_delta = 0.999,
                                                                                       max_treedepth = 15))

#Estimate anatoxins in TM mats
fit.atx.mat <-  stan(file = "HAB_toxins_Within_Mat.stan", data = model.atx.mat, chains = 3, iter = 6000,
                 warmup = 3000, refresh=100, init = init_fun_atx, control = list(adapt_delta = 0.999,
                                                            max_treedepth = 15))

 #Model checks and evaluation
library(shinystan)
library(bayesplot)
library(ggplot2)
library(rstantools)

#Can check posterior graphs in shinystan
shinystan::launch_shinystan(as.shinystan(fit.atx.mat))

#Model Checks: Within-Mat
mcmc_intervals(
  as.array(fit.atx.mat),
  pars = c("Ntheta", "Ptheta", "Atheta")) 
mcmc_intervals(
  as.array(fit.atx.mat),
  pars = c("Dtheta", "Ttheta", "Ctheta", "Rtheta")) 
mcmc_intervals(
  as.array(fit.atx.mat),
  pars = c("Beta1", "Beta2", "Beta3") )
mcmc_intervals(
  as.array(fit.atx.mat),
  pars = c("BetaAna", "BetaEpi", "BetaGeit") )
#Extract log-likelihood
log_lik_mat <- extract_log_lik(fit.atx.mat, parameter_name = "log_lik")
#Calculate WAIC
loo(log_lik_mat)
#When used waic(), "24 (58.5%) p_waic estimates greater than 0.4. We recommend trying loo instead."

#Model Checks: River-Wide
mcmc_intervals(
  as.array(fit.atx.river),
  pars = c("Ntheta", "Ptheta", "Atheta")) 
mcmc_intervals(
  as.array(fit.atx.river),
  pars = c("Dtheta", "Ttheta", "Ctheta", "Rtheta")) 
mcmc_intervals(
  as.array(fit.atx.river),
  pars = c("Beta1", "Beta2", "Beta3", "Beta4") )
mcmc_intervals(
  as.array(fit.atx.river),
  pars = c("BetaGreen", "BetaMicro", "BetaAna", "BetaNFix") )



#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.atx.river, permuted=FALSE), 
        file = here::here("data/Outputs for Obs vs Real/Anatoxin_Riverwide.rds"))
saveRDS(rstan::extract(fit.atx.mat, permuted=FALSE), 
        file = here::here("data/Outputs for Obs vs Real/Anatoxin_Withinmat.rds"))

#For building the latent state vs predictions plots
saveRDS(rstan::extract(fit.atx.river, pars = c('Beta0', 'Beta1', 'Beta2', 'Beta3', 'Beta4',
                                              'BetaGreen','BetaMicro', 'BetaAna', 
                                              'BetaNFix', 'Ntheta','Ptheta', 
                                              'Atheta', 'Dtheta', 'Ttheta', 
                                              'Ctheta', 'Rtheta', 'sigma_p',
                                              'tox_raw')), 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_River_predictions.rds"))
saveRDS(rstan::extract(fit.atx.mat, pars = c('Beta0', 'Beta1', 'Beta2', 'Beta3',
                                             'Ntheta','Ptheta', 'Atheta', 'Dtheta',
                                             'Ttheta', 'Ctheta', 'Rtheta', 'sigma_p',
                                             'tox_raw')), 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_Mat_predictions.rds"))



##### Exploratory business: Significant lag between toxins and abundances?#####
lag_df <- data.frame(time = TM_latent$time,
                     ATX = anatoxin_data$ATX_all_ug_g/1000, 
                     Anabaena_mat = exp(TM_latent$Anabaena),
                     Microcoleus_river = exp(alltaxatime$microcoleus)[-c(14:15, 29:30)],
                     Anabaena_river = exp(alltaxatime$anabaena_cylindrospermum[-c(14:15, 29:30)]))

#Plots
ccf(lag_df$Anabaena_mat, lag_df$ATX, lag.max = 5)
ccf(lag_df$Anabaena_river, lag_df$ATX, lag.max = 5)
ccf(lag_df$Microcoleus_river, lag_df$ATX, lag.max = 5)
  
                      
