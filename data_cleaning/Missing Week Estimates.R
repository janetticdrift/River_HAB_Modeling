#Missing week estimates
#This code prepares data for using with STAN to estimate the percent cover values and
#microscopy proportion values in the off-weeks of 2022 and 2024.

#Packages----
library(rstan)
library(tidyverse)
library(dataRetrieval)

#Read in needed data
source(here::here("data_cleaning/cleaning_HAB.R"))

#Tidy dataframes into format needed for STAN

#####
#This dataframe (yeardata) is created for HAB_all_years.stan, HAB_abiotic.stan, 
#and HAB_biotic.stan
#####
weekdata <- cover_indexweek %>% 
  dplyr::select(-c(timestep, field_date)) %>% 
  group_by(year) %>% 
  complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% 
  replace(is.na(.), -99) %>% 
  ungroup() %>% 
  dplyr::filter(year == "2022") %>% 
  mutate(reach = as.numeric(factor(reach))) %>% 
  #mutate(across(green_algae:other_nfixers, round, 0)) %>% #Round numbers to no decimal places
  mutate(across(everything(), ~replace(., . == 0, 1))) #Cannot have zeros for log transforming

# Eventually I should log-transform data here first, when cleaning up code
#Fills in missing weeks for years that sampled bimonthly, and sets missing entries to -99

yeardata <- cover_indexweek %>% 
  dplyr::select(-c(timestep, field_date)) %>% 
  group_by(year) %>% 
  complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% 
  replace(is.na(.), -99) %>% 
  ungroup() %>% 
  mutate(reach = as.numeric(factor(reach))) %>% 
  #mutate(across(green_algae:other_nfixers, round, 0)) %>% #Round numbers to no decimal places
  mutate(across(everything(), ~replace(., . == 0, 1))) #Cannot have zeros for log transforming

#####
#This dataframe (yearmatdata) is created for HAB_mat_community.stan
#####

yearmatdata_TM <- micro_indexweek %>% 
  dplyr::select(-c(timestep, location, field_date, slide_rep, date_analyzed, method)) %>%
  filter(sample_type == "TM") %>% 
  group_by(year) %>%
  complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% 
  mutate(sample_type = replace_na(sample_type, "TM")) %>% 
  replace(is.na(.), -99)

yearmatdata_TAC <- micro_indexweek %>% 
  dplyr::select(-c(timestep, location, field_date, slide_rep, date_analyzed, method)) %>%
  filter(sample_type == "TAC") %>% 
  group_by(year) %>%
  complete(nesting(site_reach, site, reach), week = seq(1, max(week), 1L)) %>% 
  mutate(sample_type = replace_na(sample_type, "TAC")) %>% 
  replace(is.na(.), -99)

yearmatdata <- rbind(yearmatdata_TM, yearmatdata_TAC) %>% 
  mutate(across(everything(), ~replace(., . == 0, 1))) #Cannot have zeros for log transforming

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
  mutate(firstday = if_else(week == 1 & (year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday, bare_biofilm) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) %>% 
  mutate(across(green_algae:other_nfixers, log)) %>%
  mutate(across(everything(), ~replace(.x, is.nan(.x), -99)))
#mutate(across(1:5, round, 0)) #green_algae:microcoleus to pull out, comment out if logging


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
#MULTI SPECIES and MULTI-REACH - Gather data into STAN list format
library(abind)

#Clean and transform into a 2D dataframe
mattaxareach <- yearmatdata %>% 
  mutate(firstday = if_else(week == 1 & (year == 2023), 1, 0))

for(i in 1:nrow(mattaxareach)){
  if(microcoleus == -99) {
   firstsample = 0 
  } else {
    firstsample = 1
  }
}

  filter(year == 2022) %>% #TEST ONE YEAR FIRST
  unite("uniqueID", c(year, week), sep = "_") %>% 
  relocate(firstday) %>% 
  mutate(across(anabaena_and_cylindrospermum:rare, log)) %>% #logtransform
  mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) %>%  #reset the -99s
  select(!c(site_reach, site)) %>% 
  select(c(1:5, 10)) %>%  #TEST 2 SPECIES FIRST
  filter(sample_type == "TM")  #Evaluate TM and TAC separatedly

#Convert dataframe into a list
spreach <- split(mattaxareach, mattaxareach$reach) #weeks, species, reaches


model.1 <- list("uniqueID" = nrow(spreach[["1S"]]), 
                "Nreach" = length(spreach),
                "Nspecies" = as.integer(ncol(spreach[["1S"]])-4),#take out first 4 col: firstday:sample_type
                "N" = lapply(spreach, function(df) df[, -c(1:4)]) #take out first 3 colums
)


#"firstdays" = alltaxatime$firstday,
#"id" = c(1,1,1,1),
#-------------------------------------------------------------------------------------------------
#Run models

setwd(here::here("data_cleaning")) #Set working directory to current folder

options(mc.cores = parallel::detectCores())
#####TARGET MATS
#All years, one species, 3 reaches
fit.m1 <-  stan(file = "HAB_mat_community.stan", data = model.1, chains = 3, iter = 10000,
                warmup = 5000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 20))


######RIVER WIDE

#All years, all species, averaged reach
fit.m4 <-  stan(file = "HAB_all_years.stan", data = model.4, chains = 3, iter = 10000,
                warmup = 5000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 20))

#Only biotic variables, all species, averaged reach
fit.m5 <-  stan(file = "HAB_biotic.stan", data = model.4, chains = 3, iter = 10000,
                warmup = 5000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 20))

#Only abiotic variables, all species, averaged reach
fit.m6 <-  stan(file = "HAB_abiotic.stan", data = model.4, chains = 3, iter = 10000,
                warmup = 5000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 20))

#-------------------------------------------------------------------------------------------------
#Model checks and evaluation
library(shinystan)
library(bayesplot)
library(ggplot2)
shinystan::launch_shinystan(fit.m4)
print(fit.m4, par = "Ptheta")

#Save output for cleaning and visualizing in data_analysis/model_vs_real_data.R
avg.reach.output <- as.data.frame(rstan::extract(fit.m4, permuted=FALSE))


#Save raw parameter estimates
saveRDS(
  avg.reach.output, 
  file = here::here("data/Bayes_avg_reach_output.rds")
) 
#To read: object <- readRDS(here:here("data/Bayes_avg_reach_output.rds"))
