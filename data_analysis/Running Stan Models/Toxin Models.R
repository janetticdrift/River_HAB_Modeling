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

                                 ######Isolate Microcoleus mats#####

toxinsTM <- toxindf %>% 
  dplyr::filter(sample_type == "Microcoleus") %>% 
  dplyr::mutate(field_date = replace(field_date, field_date == as.Date("2023-07-11"), 
                                     as.Date("2023-07-10"))) %>%  #Replace 2023-07-11 with 07/10 so they are on the same week
  dplyr::mutate(field_date = replace(field_date, field_date == as.Date("2022-09-06"),
                                     as.Date("2022-09-08"))) %>%  #Replace 2022/09/06 to 09/08 so they are on the same week 
  dplyr::filter(field_date != as.Date("2024-06-19")) %>% #Remove 6/19/24 as it wasn't sampled in microscopy data
  dplyr::mutate(field_date = as.Date(field_date, format="%m/%d/%y")) 

#Save cleaned output for visualizing observational vs latent states
saveRDS(toxinsTM, 
        file = here::here("data/Outputs for Obs vs Real/obs_toxins_TM.rds"))


anaCsplit <- toxinsTM %>% 
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

atx_TM <- rbind(year1_index, year2_index, year3_index) %>% 
  dplyr::select(-c(field_date, sample_type, timestep)) %>% 
  group_by(year) %>% 
  complete(nesting(reach), week = seq(min(week), max(week), 1L)) %>% 
  ungroup() %>%
  dplyr::mutate(reach = as.numeric(factor(reach))) 

#Gather data into stan list format
anatoxin_data_TM1 <- atx_TM %>% 
  group_by(year, week) %>% 
  dplyr::summarise(ATX_all_ug_afdm_g = mean(ATX_all_ug_afdm_g, na.rm = TRUE)) %>% #Average across reaches, removing reaches where no ATX was collected
  dplyr::mutate(ATX_all_ug_afdm_g = round(ATX_all_ug_afdm_g, digits = 3), #Editing data for poisson distribution
                ATX_all_ug_afdm_g = ATX_all_ug_afdm_g*1000) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) %>% 
  mutate(firstday = if_else(week == 1 & (year == 2022 | year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) 

#Save cleaned output for visualizing observational vs latent states
saveRDS(anatoxin_data_TM1, 
        file = here::here("data/Outputs for Sims and Model Fits/obs_toxins_data_TM.rds"))

#Replace -99s with zeros and add an "Is Observed?" column, because Poisson cannot take negatives
anatoxin_data_TM <- anatoxin_data_TM1 %>%  
dplyr::mutate(is_obs  = ifelse(ATX_all_ug_afdm_g == -99, 0, 1), #Editing data for poisson
              ATX_all_ug_afdm_g = ifelse(ATX_all_ug_afdm_g == -99, 0, ATX_all_ug_afdm_g))


                         ######Isolate Toxins from Anabaena mats#####

toxinsTAC <- toxindf %>% 
  dplyr::filter(sample_type == "Anabaena") %>% 
  dplyr::mutate(field_date = replace(field_date, field_date == as.Date("2022-09-06"),
                                     as.Date("2022-09-08"))) %>%  #Replace 2022/09/06 to 09/08 so they are on the same week 
  dplyr::mutate(field_date = as.Date(field_date, format="%m/%d/%y")) 

#Save cleaned output for visualizing observational vs latent states
saveRDS(toxinsTAC, 
        file = here::here("data/Outputs for Obs vs Real/obs_toxins_TAC.rds"))


anaCsplit <- toxinsTAC %>% 
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
                                 timestep == 6 ~ 11))

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
                                 timestep == 4 ~ 7))

atx_TAC <- rbind(year1_index, year2_index, year3_index) %>% 
  dplyr::select(-c(field_date, sample_type, timestep)) %>% 
  group_by(year) %>% 
  complete(nesting(reach), week = seq(min(week), max(week), 1L)) %>% 
  ungroup() %>%
  dplyr::mutate(reach = as.numeric(factor(reach))) 

#Gather data into stan list format
anatoxin_data_TAC1 <- atx_TAC %>% 
  group_by(year, week) %>% 
  dplyr::summarise(ATX_all_ug_afdm_g = mean(ATX_all_ug_afdm_g, na.rm = TRUE)) %>% #Average across reaches, removing reaches where no ATX was collected
  dplyr::mutate(ATX_all_ug_afdm_g = round(ATX_all_ug_afdm_g, digits = 3), #Editing data for poisson distribution
                ATX_all_ug_afdm_g = ATX_all_ug_afdm_g*1000) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) %>% 
  mutate(firstday = if_else(week == 1 & (year == 2022 | year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) 

#Save cleaned output for visualizing observational vs latent states
saveRDS(anatoxin_data_TAC1, 
        file = here::here("data/Outputs for Sims and Model Fits/obs_toxins_data_TAC.rds"))

#Replace -99s with zeros and add an "Is Observed?" column, because Poisson cannot take negatives
anatoxin_data_TAC <- anatoxin_data_TAC1 %>%  
  dplyr::mutate(is_obs  = ifelse(ATX_all_ug_afdm_g == -99, 0, 1), #Editing data for poisson
                ATX_all_ug_afdm_g = ifelse(ATX_all_ug_afdm_g == -99, 0, ATX_all_ug_afdm_g))


#---------------------------------------------------------------------------------------
#CREATE MODEL FOR MICROCOLEUS WITHIN-MAT MICROSCOPY DATA

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

#Create design matrix for toxins from Microcoleus mats
X1TM <- cbind(
  intercept = 1,
  TM_latent[, -1],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)],
  temp = stand_nut$temp_C[-c(14:15, 29:30)],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)],
  rad = swradiation$stand_rad[-c(14:15, 29:30)]
  # DIN = stand_nut$DIN[-c(14:15, 29:30)]
)

#Combine with other information into model list
model.atx.matTM <- list("uniqueID" = nrow(anatoxin_data_TM),
                  "is_obs" = anatoxin_data_TM$is_obs, #poisson edit
                  "firstdays" = anatoxin_data_TM$firstday,
                  "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                  "Nspecies" = as.integer(ncol(TM_latent)-1),
                  "X1" = X1TM,
                  "Npredictors" = ncol(X1TM)
)

#---------------------------------------------------------------------------------------
#Read in the model lists for each of the toxin model variations (All Variables, Biotic,
  #Abiotic, Abiotic No Nutrients)
source(here::here("data_analysis/Running Stan Models/Toxin Models_Cleaning Latent Algal Abundances.R"))

  #Objects read in are in the format of "model.atx.river[Target Mat]_[Model Name Here]

#Use toxins samples from Anabaena mats
model.atx.riverTAC <- list("uniqueID" = nrow(anatoxin_data_TAC),
                        "is_obs" = anatoxin_data_TAC$is_obs, #poisson edit, was this an observed day?
                        "firstdays" = anatoxin_data_TAC$firstday,
                        "Toxins" = as.integer(anatoxin_data_TAC$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                        "Nspecies" = as.integer(ncol(River_latent)-1),
                        "X2" = X2TAC,
                        "Npredictors" = ncol(X2TAC)
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
  Beta3 = 0)
#tox_nc = rep(0, nrow(anatoxin_data_TM)) #nrow is the time length
 

#Estimate anatoxins using river-wide assemblages
    #Toxins from Microcoleus mats
#All Variables
fit.atx.riverTMall <-  stan(file = "HAB_toxins_River_Wide.stan", data = model.atx.riverTM_All, chains = 3, iter = 6000,
                       warmup = 3000, refresh=100, init = init_fun_atx, control = list(adapt_delta = 0.999,
                                                                                       max_treedepth = 15))

#Biotic Only
fit.atx.riverTMbiotic <-  stan(file = "HAB_toxins_River_Wide.stan", data = model.atx.riverTM_Biotic, chains = 3, iter = 6000,
                         warmup = 3000, refresh=100, init = init_fun_atx, control = list(adapt_delta = 0.999,
                                                                                         max_treedepth = 15))

#Abiotic Only
fit.atx.riverTMabiotic <-  stan(file = "HAB_toxins_River_Wide.stan", data = model.atx.riverTM_Abiotic, chains = 3, iter = 6000,
                               warmup = 3000, refresh=100, init = init_fun_atx, control = list(adapt_delta = 0.999,
                                                                                               max_treedepth = 15))

#Abiotic Minus Nutrients
fit.atx.riverTMabioticnonut <-  stan(file = "HAB_toxins_River_Wide.stan", data = model.atx.riverTM_AbioticNoNut, chains = 3, iter = 6000,
                               warmup = 3000, refresh=100, init = init_fun_atx, control = list(adapt_delta = 0.999,
                                                                                               max_treedepth = 15))


     #Toxins from Anabaena mats
fit.atx.riverTAC <-  stan(file = "HAB_toxins_River_Wide.stan", data = model.atx.riverTAC, chains = 3, iter = 6000,
                       warmup = 3000, refresh=100, init = init_fun_atx, control = list(adapt_delta = 0.999,
                                                                                       max_treedepth = 15))

#Estimate anatoxins using TM microscopy assemblages
fit.atx.mat <-  stan(file = "HAB_toxins_Within_Mat.stan", data = model.atx.matTM, chains = 3, iter = 6000,
                 warmup = 3000, refresh=100, init = init_fun_atx, control = list(adapt_delta = 0.999,
                                                            max_treedepth = 15))



#For building the observation vs latent state plots
                              ###River-Wide###
  #Microcoleus Mat Anatoxins
saveRDS(rstan::extract(fit.atx.riverTMall, permuted=FALSE), 
        file = here::here("data/Outputs for Obs vs Real/Anatoxin_TM_Riverwide_All.rds"))
saveRDS(rstan::extract(fit.atx.riverTMbiotic, permuted=FALSE), 
        file = here::here("data/Outputs for Obs vs Real/Anatoxin_TM_Riverwide_Biotic.rds"))
saveRDS(rstan::extract(fit.atx.riverTMabiotic, permuted=FALSE), 
        file = here::here("data/Outputs for Obs vs Real/Anatoxin_TM_Riverwide_Abiotic.rds"))
saveRDS(rstan::extract(fit.atx.riverTMabioticnonut, permuted=FALSE), 
        file = here::here("data/Outputs for Obs vs Real/Anatoxin_TM_Riverwide_AbioticNoNut.rds"))

  #Anabaena Mat Anatoxins
saveRDS(rstan::extract(fit.atx.riverTAC, permuted=FALSE), 
        file = here::here("data/Outputs for Obs vs Real/Anatoxin_TAC_Riverwide.rds"))

                              ###Within-Mat###
saveRDS(rstan::extract(fit.atx.mat, permuted=FALSE), 
        file = here::here("data/Outputs for Obs vs Real/Anatoxin_Withinmat.rds"))

#For building the latent state vs predictions plots
                              ###River-Wide###
saveRDS(rstan::extract(fit.atx.riverTMall, pars = c('Beta0', 'Beta1', 'Beta2', 'Beta3', 'Beta4',
                                              'Phi0',
                                              'PhiAna', 
                                              'Ntheta',
                                              'Ptheta',
                                              'Atheta',
                                              'Dtheta', 'Ttheta', 
                                              'Ctheta', 'Rtheta', 'sigma_p',
                                              'tox_raw')), 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_All_predictions.rds"))
saveRDS(rstan::extract(fit.atx.riverTMbiotic, pars = c('Beta0', 'Beta1', 'Beta2', 'Beta3', 'Beta4',
                                                    'Phi0',
                                                    'PhiAna', 
                                                    'Ntheta',
                                                    'Ptheta',
                                                    'Atheta',
                                                    'Dtheta', 'Ttheta', 
                                                    'Ctheta', 'Rtheta', 'sigma_p',
                                                    'tox_raw')), 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_Biotic_predictions.rds"))
saveRDS(rstan::extract(fit.atx.riverTMabiotic, pars = c('Beta0', 'Beta1', 'Beta2', 'Beta3', 'Beta4',
                                                       'Phi0',
                                                       'PhiAna', 
                                                       'Ntheta',
                                                       'Ptheta',
                                                       'Atheta',
                                                       'Dtheta', 'Ttheta', 
                                                       'Ctheta', 'Rtheta', 'sigma_p',
                                                       'tox_raw')), 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_Abiotic_predictions.rds"))
saveRDS(rstan::extract(fit.atx.riverTMabioticnonut, pars = c('Beta0', 'Beta1', 'Beta2', 'Beta3', 'Beta4',
                                                       'Phi0',
                                                       'PhiAna', 
                                                       'Ntheta',
                                                       'Ptheta',
                                                       'Atheta',
                                                       'Dtheta', 'Ttheta', 
                                                       'Ctheta', 'Rtheta', 'sigma_p',
                                                       'tox_raw')), 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TM_River_AbioticNoNut_predictions.rds"))



saveRDS(rstan::extract(fit.atx.riverTAC, pars = c('Beta0', 'Beta1', 'Beta2', 'Beta3', 'Beta4',
                                                 'Phi0',
                                                 'PhiAna', 
                                                 'Ntheta',
                                                 'Ptheta',
                                                 'Atheta',
                                                 'Dtheta', 'Ttheta', 
                                                 'Ctheta', 'Rtheta', 'sigma_p',
                                                 'tox_raw')), 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_TAC_River_predictions.rds"))

                                      ###Within-Mat###
saveRDS(rstan::extract(fit.atx.mat, pars = c('Beta0', 'Beta1', 'Beta2', 'Beta3',
                                             'Ntheta',
                                             'Ptheta',
                                             'Atheta',
                                             'Dtheta', 'Ttheta', 
                                             'Ctheta', 'Rtheta', 'sigma_p', 'phi',
                                             'tox_raw')), 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Anatoxin_Mat_predictions.rds"))


#  #Model checks and evaluation
# library(shinystan)
# library(bayesplot)
# library(ggplot2)
# library(rstantools)
# 
# #Can check posterior graphs in shinystan
shinystan::launch_shinystan(as.shinystan(fit.atx.riverTMall))
# 
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
  pars = c("phi") )
#Extract log-likelihood
log_lik_mat <- extract_log_lik(fit.atx.riverTMall, parameter_name = "log_lik")
#Calculate WAIC
loo(log_lik_mat)
waic(log_lik_mat)
#When used waic(), "24 (58.5%) p_waic estimates greater than 0.4. We recommend trying loo instead."

# #Model Checks: River-Wide
mcmc_intervals(
  as.array(fit.atx.riverTMall),
  pars = c("Ntheta", "Ptheta", "Atheta"))
mcmc_intervals(
  as.array(fit.atx.riverTMall),
  pars = c("Dtheta", "Ttheta", "Ctheta", "Rtheta"))
mcmc_intervals(
  as.array(fit.atx.riverTMall),
  pars = c("Beta1", "Beta2", "Beta3", "Beta4") )
mcmc_intervals(
  as.array(fit.atx.riverTMall),
  pars = c("PhiAna") )

##### Exploratory business: Significant lag between toxins and abundances? Or Env drivers?#####
#Create new dataframe of variables that we want to check the lagged relationship between
lag_df <- data.frame(time = TM_latent$time,
                     ATX = anatoxin_data_TM$ATX_all_ug_afdm_g/1000, 
                     Anabaena_mat = exp(TM_latent$Anabaena),
                     Microcoleus_river = exp(alltaxatime$microcoleus)[-c(14:15, 29:30)],
                     Anabaena_river = exp(alltaxatime$anabaena_cylindrospermum)[-c(14:15, 29:30)],
                     nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)],
                     phosphate = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)],
                     ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)],
                     discharge = discharge$stand_discharge[-c(14:15, 29:30)],
                     conduc = stand_nut$cond_uS_cm[-c(14:15, 29:30)],
                     temp = stand_nut$temp_C[-c(14:15, 29:30)],
                     radiaion = swradiation$stand_rad[-c(14:15, 29:30)]
                     )

#Cross-Correlation
ccf(lag_df$Anabaena_mat, lag_df$ATX, lag.max = 5)
ccf(lag_df$Anabaena_river, lag_df$ATX, lag.max = 5)
ccf(lag_df$Microcoleus_river, lag_df$ATX, lag.max = 5)
 ccf(lag_df$nitrate, lag_df$ATX, lag.max = 5)

#Instead of CCF plots, stack regression relationship 
#Create dataframe of different time lags
lag_df_regression <- lag_df %>%
  dplyr::mutate("Anabaena (t-1)" = lag(Anabaena_river, 1),     #This column lags the abundance of Ana by 1 time step
                "Anabaena (t-2)" = lag(Anabaena_river, 2),     #This column lags the abundance of Anabaena by 2 time steps
                "Anabaena (t-3)" = lag(Anabaena_river, 3),     #and so on...
                "Microcoleus (t-1)" = lag(Microcoleus_river, 1),     
                "Microcoleus (t-2)" = lag(Microcoleus_river, 2),
                "Microcoleus (t-3)" = lag(Microcoleus_river, 3)) %>% 
  dplyr::select(!Anabaena_mat) %>% #Remove within-mat microscopy abundance
  dplyr::rename("Anabaena (t)" = Anabaena_river, "Microcoleus (t)" = Microcoleus_river)

#Pivot dataframe to group by amount of lag
lag_dfpivot <- lag_df_regression %>% 
  pivot_longer(cols = "Microcoleus (t)":"Microcoleus (t-3)",
               names_to = "lag",
               values_to = "Abundance") %>% 
  #Add a column for categorizing the Microcoleus from Anabaena samples
  dplyr::mutate(taxa = case_when(
                 lag = str_starts(lag, "Micro") ~ "Microcoleus",
                 lag = str_starts(lag, "Ana") ~ "Anabaena"))

#Plot lagged relationships between percent cover abundances of Microcoleus and Anabaena
  #and toxin concentrations

microlag <- ggplot(subset(lag_dfpivot, taxa %in% "Microcoleus"), aes(x = Abundance, y = ATX, color = lag)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("black", "#F1C6C6", "#DE7C7C", "#A52A29"),
                     breaks = c("Microcoleus (t)",
                                "Microcoleus (t-1)",
                                "Microcoleus (t-2)",
                                "Microcoleus (t-3)")) +
  coord_cartesian(y = c(0, 46)) + 
  labs(x = "Percent Cover Abundance", y = "Toxin Concentration (ug/g)", color = "Lag Time") +
  theme_bw()

analag <- ggplot(subset(lag_dfpivot, taxa %in% "Anabaena"), aes(x = Abundance, y = ATX, color = lag)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("black", "#B7E5E5", "#66a4a2", "#006664"),
                     breaks = c("Anabaena (t)",
                                "Anabaena (t-1)",
                                "Anabaena (t-2)",
                                "Anabaena (t-3)")) +
  coord_cartesian(y = c(0, 46)) + 
  labs(x = "Percent Cover Abundance", y = "Toxin Concentration (ug/g)", color = "") +
  theme_bw()

(microlag / analag) +
  plot_layout(guides = "collect", axes = "collect") + 
  plot_annotation(title = "Toxins from TM Mats")

#Calculate correlation between toxin concentration and Anabanena at t-2
cor(lag_df_regression$ATX, lag_df_regression$`Anabaena (t-2)`, use="complete.obs")
