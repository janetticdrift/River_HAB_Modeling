#Missing week estimates
#This code prepares data for using with STAN to estimate the percent cover values and
#microscopy proportion values in the off-weeks of 2022 and 2024.

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
#This dataframe (yeardata) is created for HAB_all_years.stan, HAB_abiotic.stan, 
#and HAB_biotic.stan
#####
# weekdata <- cover_indexweek %>% 
#   dplyr::select(-c(timestep, field_date)) %>% 
#   group_by(year) %>% 
#   dplyr::mutate(across(green_algae:other_nfixers,
#                        ~ . + pseudocount)) %>% #Cannot have zeros for log transforming
#   complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% 
#   replace(is.na(.), -99) %>% 
#   ungroup() %>% 
#   #dplyr::filter(year == "2022") %>% 
#   dplyr::mutate(reach = as.numeric(factor(reach))) 
#   #mutate(across(green_algae:other_nfixers, round, 0)) %>% #Round numbers to no decimal places

# Eventually I should log-transform data here first, when cleaning up code
#Fills in missing weeks for years that sampled bimonthly, and sets missing entries to -99

yeardata <- cover_indexweek %>% 
  dplyr::select(-c(timestep, field_date)) %>% 
  group_by(year) %>% 
  dplyr::mutate(across(green_algae:other_nfixers,
                       ~ . + pseudocount)) %>% #Cannot have zeros for log transforming
  complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% 
  replace(is.na(.), -99) %>% 
  ungroup() %>% 
  dplyr::mutate(reach = as.numeric(factor(reach))) 
  #mutate(across(green_algae:other_nfixers, round, 0)) %>% #Round numbers to no decimal places

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

yearmatdata <- rbind(yearmatdata_TM, yearmatdata_TAC) %>% 
  arrange(year, week)

#-------------------------------------------------------------------------------------------------
# #SINGLE SPECIES - Gather data into STAN list format
# 
# #This dataframe uses HAB_two_species.stan!
# 
# #Change formatting of only green algae to start
# green_algae <- weekdata %>% 
#   dplyr::select(reach, green_algae) %>% 
#   mutate(row = rep(seq(1:13), length(unique(reach)))) %>% #13 = number of collection days
#   pivot_wider(names_from = reach, values_from = green_algae) %>% 
#   select(-row)
# 
# #Finished dataframe: row is # of weeks, columns is reach number
# 
# model.1 <- list("Nweeks" = length(unique(weekdata$week)), 
#                 "Nreach" = length(unique(weekdata$site_reach)),
#                 "N" = green_algae #Only one species right now
# )

#-------------------------------------------------------------------------------------------------
# #TWO SPECIES - Gather data into STAN list format
# 
# #Currently includes all taxa in this code chunk
# green_micro <- weekdata %>% 
#   group_by(week) %>% 
#   dplyr::summarise(green_algae = mean(green_algae), microcoleus = mean(microcoleus),
#                    anabaena_cylindrospermum = mean(anabaena_cylindrospermum),
#                    other_nfixers = mean(other_nfixers)) %>% #Bare is not a living species
#   mutate_if(is.numeric, log) %>%
#   mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) %>%
#   select(-week) 
# #mutate(across(1:5, round, 0)) #green_algae:microcoleus to pull out, comment out if logging
# 
# 
# model.2 <- list("Nweeks" = nrow(green_micro), 
#                 "Nspecies" = ncol(green_micro),
#                 "N" = green_micro #green algae and microcoleus
# )

#-------------------------------------------------------------------------------------------------
# #TWO SPECIES and MULTI-REACH - Gather data into STAN list format
# library(abind)
# 
# #Clean and transform in a 2D dataframe
# temp.spreach <- weekdata %>% 
#   select(-c(1:3, 5, 8:10)) %>% 
#   mutate(across(green_algae:microcoleus, round, 0)) 
# # mutate(across(.cols = c("green_algae":3), .fns = log)) %>%  #logtransform
# # mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) #reset the -99s
# 
# #Split data into an array by reach, then drop the reach column
# spreach.array <- abind(split(temp.spreach[, -1], temp.spreach$reach), along = 3)
# 
# #Convert array into a list
# spreach = plyr::alply(spreach.array,3, .dims = TRUE)
# 
# 
# model.3 <- list("Nweeks" = nrow(spreach[["1"]]), 
#                 "Nreach" = length(spreach),
#                 "Nspecies" = ncol(spreach[["1"]]),
#                 "N" = spreach
# )

#-------------------------------------------------------------------------------------------------
#MULTISPECIES and MULTI-YEAR - Gather data into STAN list format

alltaxatime <- yeardata %>% 
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


model.4 <- list("uniqueID" = nrow(alltaxatime), 
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
#MAT COMMUNITY PER REACH-----MICROCOLEUS
#Split data into reach subsets, and run model separately per reach

#Important note: Only Ana, Epithemia, and Geit are retained here. Their abundances
#have been re-relativized to sum to 100%

#Target Microcoleus, averaged reaches
matalltaxaM <- yearmatdata %>% 
  dplyr::filter(sample_type == "TM") %>% 
  dplyr::select(c(1:9)) %>% #Retain only Ana, Epi, and Geit
  # rowwise() %>% #Re-relativize one row at a time
  # dplyr::mutate(across(Anabaena:Geitlerinema, ~ .x / sum(c_across(Anabaena:Geitlerinema)) * 100)) %>% #Divide the abundances by the new row total, *100
  # ungroup() %>% 
  dplyr::mutate(across(c(Anabaena:Geitlerinema),
                       ~ . + pseudocount)) %>%  #Cannot have zeros for log transforming
  group_by(year, week) %>%
  dplyr::summarise(across(c(Anabaena:Geitlerinema), mean, na.rm = TRUE)) %>% 
  mutate(firstday = if_else(week == 1 & (year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) %>% 
  dplyr::mutate(across(Anabaena:Geitlerinema, log)) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), -99)))
  # dplyr::rename(Anabaena = anabaena_and_cylindrospermum, 
  #               'Epithemia Diatoms' = e_diatoms,
  #               Geitlerinema = geitlerinema,
  #               'Green Algae' = green_algae, 
  #               Microcoleus = microcoleus,
  #               'Non-Epithemia Diatoms' = non_e_diatoms,
  #               Nostoc = nostoc,
  #               'Other Coccoids' = other_coccoids,
  #               Rare = rare)


model.1 <- list("uniqueID" = nrow(matalltaxaM),
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
#MAT COMMUNITY PER REACH-----ANABAENA
#Split data into reach subsets, and run model separately per reach

#Target Anabaena, averaged reaches
matalltaxaA <- yearmatdata %>% 
  dplyr::filter(sample_type == "TAC") %>% 
  group_by(year, week) %>%
  dplyr::summarise(across(c(anabaena_and_cylindrospermum:rare), mean, na.rm = TRUE)) %>% 
  mutate(firstday = if_else(week == 2 & year == 2023 | week == 3 & year == 2024, 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) %>% 
  dplyr::mutate(across(anabaena_and_cylindrospermum:rare, log)) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) %>% 
  dplyr::rename(Anabaena = anabaena_and_cylindrospermum, 
                'Epithemia Diatoms' = e_diatoms,
                Geitlerinema = geitlerinema,
                'Green Algae' = green_algae, 
                Microcoleus = microcoleus,
                'Non-Epithemia Diatoms' = non_e_diatoms,
                Nostoc = nostoc,
                'Other Coccoids' = other_coccoids,
                Rare = rare) %>% 


model.2 <- list("uniqueID" = nrow(matalltaxaA),
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


#####TARGET MATS
#All years, one species, 3 reaches

init_fun_M <- function() list(
  sigma_p = rep(0.5, 3),     #9 is number of species in mat datasets
  sigma_o = rep(0.5, 3),
  Alpha   = rep(0, 3),
  Beta_diag = rep(0, 0, 3),     # small start
  Beta_off = matrix(0, 3, 3),
  n_nc = matrix(0, 3, nrow(matalltaxaM))
)

init_fun_A <- function() list(
  sigma_p = rep(0.5, 9),     #9 is number of species in mat datasets
  sigma_o = rep(0.5, 9),
  Alpha   = rep(0, 9),
  Beta_diag = rep(0, 0, 9),     # small start
  Beta_off = matrix(0, 9, 9),
  n_nc = matrix(0, 9, nrow(matalltaxaA))
)

#Averaged, TM
fit.m1 <-  stan(file = "HAB_mat_community.stan", data = model.1, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, init = init_fun_M, control = list(adapt_delta = 0.999,
                                                           max_treedepth = 15))
#Averaged, TAC
fit.m2 <-  stan(file = "HAB_mat_community.stan", data = model.2, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, init = init_fun_A, control = list(adapt_delta = 0.999,
                                                                            max_treedepth = 15))

######RIVER WIDE

#All years, all species, averaged reach
fit.m4 <- stan(file = "HAB_all_years.stan", data = model.4, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 13))

#Only biotic variables, all species, averaged reach
fit.m5 <-  stan(file = "HAB_biotic.stan", data = model.4, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 13))

#Only abiotic variables, all species, averaged reach
fit.m6 <-  stan(file = "HAB_abiotic.stan", data = model.4, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 13))

#Only abiotic variables, no nutrients, all species, averaged reach
fit.m7 <-  stan(file = "HAB_abiotic_nonut.stan", data = model.4, chains = 3, iter = 6000,
                warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 13))

 #-------------------------------------------------------------------------------------------------
#Model checks and evaluation
library(shinystan)
library(bayesplot)
library(ggplot2)
library(rstantools)

#Can check posterior graphs in shinystan
shinystan::launch_shinystan(fit.m1)
print(fit.m4, par = "Ptheta")


#Save within-mat output for cleaning and visualizing in data_analysis/model_vs_real_data.R

saveRDS(rstan::extract(fit.m1, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta'),
                       permuted=FALSE), 
        file = here::here("data/WithinMat_Micro.rds"))
saveRDS(rstan::extract(fit.m1, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/WithinMat_Micro_predictions.rds"))


saveRDS(rstan::extract(fit.m2, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta'),
                       permuted=FALSE), 
        file = here::here("data/WithinMat_Ana.rds"))
saveRDS(rstan::extract(fit.m2, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/WithinMat_Ana_predictions.rds"))


################
#Save river-wide output for cleaning and visualizing in data_analysis/model_vs_real_data.R

#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.m4, permuted=FALSE), 
        file = here::here("data/Riverwide_AllVariables.rds"))
#For building the latent state vs predictions plots
saveRDS(rstan::extract(fit.m4, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/Riverwide_AllVar_predictions.rds"))

#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.m5, permuted=FALSE), 
        file = here::here("data/Riverwide_Biotic.rds"))
#For building the latent state vs predictions plots, and calculating fit indices
saveRDS(rstan::extract(fit.m5, pars = c('Alpha', 'Beta', 'n', 'sigma_p')), 
        file = here::here("data/Riverwide_Biotic_predictions.rds"))


#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.m6, permuted=FALSE), 
        file = here::here("data/Riverwide_Abiotic.rds"))
#For building the latent state vs predictions plots, and calculating fit indices
saveRDS(rstan::extract(fit.m6, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Ntheta','Ptheta', 'Atheta', 
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/Riverwide_Abiotic_predictions.rds"))

#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.m7, permuted=FALSE), 
        file = here::here("data/Riverwide_Abiotic_nonut.rds"))
#For building the latent state vs predictions plots, and calculating fit indices
saveRDS(rstan::extract(fit.m7, pars = c('Alpha', 'Beta', 'n', 'sigma_p',
                                        'Dtheta', 'Ttheta', 'Ctheta', 
                                        'Rtheta')), 
        file = here::here("data/Riverwide_AbioticNonut_predictions.rds"))


#To read: object <- readRDS(here::here("data/Riverwide_AllVariables.rds"))
