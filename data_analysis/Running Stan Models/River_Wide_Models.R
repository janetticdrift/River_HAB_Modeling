###########################
#Running River-Wide Models: Estimating Latent States of Observed and Skipped Weeks
###########################
#Create dataframes used for running Stan models that estimate the latent states of each
  #algal species' percent cover abundance per week. The model also estimates the 
  #percent covers during weeks where no field observations were made in 2022 and 2024.

#RDS files saved are of model outputs, and are further used to build prediction
  #simulations, as well as calculating goodness-of-fit indices.


#Packages----
library(rstan)
library(tidyverse)
library(dataRetrieval)

#Read in needed data
source(here::here("data_cleaning/cleaning_HAB.R"))
#Set pseudocount here for tidying data before log-transforming
pseudocount <- 1

#Tidy dataframes into format needed for STAN

#####
#This dataframe (yearriverdata) is created for HAB_all_years.stan, HAB_abiotic.stan, 
#HAB_biotic.stan, and HAB_abiotic_nonut.stan
#####
yearriverdata <- cover_indexweek %>% 
  dplyr::select(-c(timestep, field_date)) %>% 
  group_by(year) %>% 
  dplyr::mutate(across(green_algae:other_nfixers,
                       ~ . + pseudocount)) %>% #Cannot have zeros for log transforming
  complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% 
  replace(is.na(.), -99) %>% 
  ungroup() %>% 
  dplyr::mutate(reach = as.numeric(factor(reach))) 


#####
#This dataframe (yearmatdata) is created for HAB_mat_community.stan
#####

yearmatdata_TM <- micro_indexweek %>% 
  dplyr::select(-c(location, field_date, slide_rep, date_analyzed, method)) %>%
  dplyr::filter(sample_type == "TM") %>% 
  pivot_wider(names_from = Species, values_from = Abundance) %>% 
  group_by(year) %>%
  complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% 
  dplyr::mutate(sample_type = replace_na(sample_type, "TM"))

yearmatdata_TAC <- micro_indexweek %>% 
  dplyr::select(-c(location, field_date, slide_rep, date_analyzed, method)) %>%
  dplyr::filter(sample_type == "TAC") %>% 
  pivot_wider(names_from = Species, values_from = Abundance) %>% 
  group_by(year) %>%
  complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% 
  dplyr::mutate(sample_type = replace_na(sample_type, "TAC"))

#Bind together target samples into one dataframe again
yearmatdata <- rbind(yearmatdata_TM, yearmatdata_TAC) %>% 
  arrange(year, week)

#-------------------------------------------------------------------------------------------------
#River-Wide - Gather percent cover data into STAN list format

alltaxatime <- yearriverdata %>% 
  group_by(year, week) %>% 
  dplyr::summarise(green_algae = mean(green_algae), microcoleus = mean(microcoleus),
                   anabaena_cylindrospermum = mean(anabaena_cylindrospermum), 
                   bare_biofilm = mean(bare_biofilm),
                   other_nfixers = mean(other_nfixers)) %>% #Average across reaches
  dplyr::mutate(firstday = if_else(week == 1 & (year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday, bare_biofilm) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) %>% 
  dplyr::mutate(across(green_algae:other_nfixers, log)) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), -99)))

#Save observed data for model fit comparisons
saveRDS(alltaxatime, 
        file = here::here("data/Model Fits/obs_river_data.rds"))

model.1 <- list("uniqueID" = nrow(alltaxatime), 
                "Nspecies" = as.integer(ncol(alltaxatime)-3),
                "firstdays" = alltaxatime$firstday,
                "id" = c(1,1,1,1),
                "N" = alltaxatime[,-(1:3)], #all species
                "nitrate" = stand_nut$nitrate_mg_N_L,
                "phos" = stand_nut$oPhos_ug_P_L,
                "ammonium" = stand_nut$ammonium_mg_N_L,
                "discharge" = discharge$stand_discharge,
                "temp" = stand_nut$temp_C,
                "cond" = stand_nut$cond_uS_cm,
                "rad" = swradiation$stand_rad
)
#-------------------------------------------------------------------------------------------------
#Within-Mat: Microcoleus - Gather percent cover data into STAN list format

#Important note: Only Ana, Epithemia, and Geit are retained here

#Target Microcoleus, averaged reaches
matalltaxaM <- yearmatdata %>% 
  dplyr::filter(sample_type == "TM") %>% 
  dplyr::select(c(1:9)) %>% #Retain only Ana, Epi, and Geit
  dplyr::mutate(across(c(Anabaena:Geitlerinema),
                       ~ . + pseudocount)) %>%  #Cannot have zeros for log transforming
  group_by(year, week) %>%
  dplyr::summarise(across(c(Anabaena:Geitlerinema), mean, na.rm = TRUE)) %>% 
  mutate(firstday = if_else(week == 1 & (year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) %>% 
  dplyr::mutate(across(Anabaena:Geitlerinema, log)) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), -99)))

#Save observed data for model fit comparisons
saveRDS(matalltaxaM, 
        file = here::here("data/Model Fits/obs_TM_data.rds"))


model.2 <- list("uniqueID" = nrow(matalltaxaM),
                "Nspecies" = as.integer(ncol(matalltaxaM)-2),#take out first 2 col: firstday and uniqueID
                "firstdays" = matalltaxaM$firstday,
                "N" = matalltaxaM[,-(1:2)],
                "nitrate" = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)], #Can subset 2024 out with 29:45
                "phos" = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)], #and also first two weeks of 2023 and 2024
                "ammonium" = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)], #Which is 14:15 and 29:30
                "discharge" = discharge$stand_discharge[-c(14:15, 29:30)],
                "temp" = stand_nut$temp_C[-c(14:15, 29:30)],
                "cond" = stand_nut$cond_uS_cm[-c(14:15, 29:30)],
                "rad" = swradiation$stand_rad[-c(14:15, 29:30)]
)


#-------------------------------------------------------------------------------------------------
#Within-Mat: Anabaena - Gather percent cover data into STAN list format

#Target Anabaena, averaged reaches
matalltaxaA <- yearmatdata %>% 
  dplyr::filter(sample_type == "TAC") %>% 
  group_by(year, week) %>%
  dplyr::summarise(across(c(Anabaena:Rare), mean, na.rm = TRUE)) %>% 
  mutate(firstday = if_else(week == 2 & year == 2023 | week == 3 & year == 2024, 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) %>% 
  dplyr::mutate(across(c(Anabaena:Rare),
                       ~ . + pseudocount)) %>%  #Cannot have zeros for log transforming
  dplyr::mutate(across(Anabaena:Rare, log)) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), -99)))


model.3 <- list("uniqueID" = nrow(matalltaxaA),
                "Nspecies" = as.integer(ncol(matalltaxaA)-2),#take out first 2 col: firstday and uniqueID
                "firstdays" = matalltaxaA$firstday,
                "N" = matalltaxaA[,-(1:2)],
                "nitrate" = stand_nut$nitrate_mg_N_L[-c(1:2, 14:16, 29:31, 39:45)], #Remove first weeks where TAC was not sampled
                "phos" = stand_nut$oPhos_ug_P_L[-c(1:2, 14:16, 29:31, 39:45)], 
                "ammonium" = stand_nut$ammonium_mg_N_L[-c(1:2, 14:16, 29:31, 39:45)],
                "discharge" = discharge$stand_discharge[-c(1:2, 14:16, 29:31, 39:45)],
                "temp" = stand_nut$temp_C[-c(1:2, 14:16, 29:31, 39:45)],
                "cond" = stand_nut$cond_uS_cm[-c(1:2, 14:16, 29:31, 39:45)],
                "rad" = swradiation$stand_rad[-c(1:2, 14:16, 29:31, 39:45)]
)

#-------------------------------------------------------------------------------------------------
#Run models

setwd(here::here("data_cleaning")) #Set working directory to current folder

options(mc.cores = parallel::detectCores())

###### River-Wide

#All years, all species, averaged reach
fit.m1.1 <- stan(file = "HAB_all_years.stan", data = model.1, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 13))

#Only biotic variables, all species, averaged reach
fit.m1.2 <-  stan(file = "HAB_biotic.stan", data = model.1, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 13))

#Only abiotic variables, all species, averaged reach
fit.m1.3 <-  stan(file = "HAB_abiotic.stan", data = model.1, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 13))

#Only abiotic variables, no nutrients, all species, averaged reach
fit.m1.4 <-  stan(file = "HAB_abiotic_nonut.stan", data = model.1, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 13))
###### Within-Mat
# Set values to initialize models

init_fun_M <- function() list(
  sigma_p = rep(0.5, 3),     #3 is number of species in TM mat datasets
  sigma_o = rep(0.5, 3),
  Alpha   = rep(0, 3),
  Beta_diag = rep(0, 0, 3),     # small start
  Beta_off = matrix(0, 3, 3),
  n_nc = matrix(0, 3, nrow(matalltaxaM))
)
init_fun_A <- function() list(
  sigma_p = rep(0.5, 9),     #9 is number of species in TAC mat datasets
  sigma_o = rep(0.5, 9),
  Alpha   = rep(0, 9),
  Beta_diag = rep(0, 0, 9),     # small start
  Beta_off = matrix(0, 9, 9),
  n_nc = matrix(0, 9, nrow(matalltaxaA))
)

#Averaged, TM
fit.m2 <-  stan(file = "HAB_mat_community.stan", data = model.2, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, init = init_fun_M, control = list(adapt_delta = 0.999,
                                                                              max_treedepth = 15))
#Averaged, TAC
fit.m3 <-  stan(file = "HAB_mat_community.stan", data = model.3, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, init = init_fun_A, control = list(adapt_delta = 0.999,
                                                                              max_treedepth = 15))

 #-------------------------------------------------------------------------------------------------
#Model checks and evaluation
library(shinystan)
library(bayesplot)
library(ggplot2)
library(rstantools)

#Can check posterior graphs in shinystan
shinystan::launch_shinystan(fit.m1)
print(fit.m4, par = "Ptheta")


################
#Save River-Wide output for cleaning and visualizing in data_analysis/model_vs_real_data.R

#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.m1.1, permuted=FALSE), 
        file = here::here("data/Riverwide_AllVariables.rds"))
#For building the latent state vs predictions plots
saveRDS(rstan::extract(fit.m1.1, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/Model Fits/Riverwide_AllVar_predictions.rds"))

#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.m1.2, permuted=FALSE), 
        file = here::here("data/Riverwide_Biotic.rds"))
#For building the latent state vs predictions plots, and calculating fit indices
saveRDS(rstan::extract(fit.m1.2, pars = c('Alpha', 'Beta', 'n', 'sigma_p')), 
        file = here::here("data/Model Fits/Riverwide_Biotic_predictions.rds"))


#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.m1.3, permuted=FALSE), 
        file = here::here("data/Riverwide_Abiotic.rds"))
#For building the latent state vs predictions plots, and calculating fit indices
saveRDS(rstan::extract(fit.m1.3, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/Model Fits/Riverwide_Abiotic_predictions.rds"))

#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.m1.4, permuted=FALSE), 
        file = here::here("data/Riverwide_Abiotic_nonut.rds"))
#For building the latent state vs predictions plots, and calculating fit indices
saveRDS(rstan::extract(fit.m1.4, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/Model Fits/Riverwide_AbioticNonut_predictions.rds"))


##Save Within-Mat output for cleaning and visualizing in data_analysis/model_vs_real_data.R

#TM output
saveRDS(rstan::extract(fit.m2, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta'), permuted=FALSE), 
        file = here::here("data/WithinMat_Micro.rds"))
saveRDS(rstan::extract(fit.m2, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/Model Fits/WithinMat_Micro_predictions.rds"))

#TAC output
saveRDS(rstan::extract(fit.m3, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta'), permuted=FALSE), 
        file = here::here("data/WithinMat_Ana.rds"))
saveRDS(rstan::extract(fit.m3, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/Model Fits/WithinMat_Ana_predictions.rds"))


#To read: object <- readRDS(here::here("data/file_name.rds"))
