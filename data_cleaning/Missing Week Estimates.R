#Missing week estimates

#Packages----
library(rstan)
library(tidyverse)
library(dataRetrieval)

#Read in needed data
source(here::here("data_cleaning/cleaning_HAB.R"))

#Tidy dataframe into format needed for STAN
#This dataframe uses HAB_missing_weeks.stan

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
#-------------------------------------------------------------------------------------------------
#SINGLE SPECIES - Gather data into STAN list format

#This dataframe uses HAB_two_species.stan!

#Change formatting of only green algae to start
green_algae <- weekdata %>% 
  dplyr::select(reach, green_algae) %>% 
  mutate(row = rep(seq(1:13), length(unique(reach)))) %>% #13 = number of collection days
  pivot_wider(names_from = reach, values_from = green_algae) %>% 
  select(-row)

#Finished dataframe: row is # of weeks, columns is reach number

model.1 <- list("Nweeks" = length(unique(weekdata$week)), 
                "Nreach" = length(unique(weekdata$site_reach)),
                "N" = green_algae #Only one species right now
)

#-------------------------------------------------------------------------------------------------
#TWO SPECIES - Gather data into STAN list format

#Currently includes all taxa in this code chunk
green_micro <- weekdata %>% 
  group_by(week) %>% 
  dplyr::summarise(green_algae = mean(green_algae), microcoleus = mean(microcoleus),
                   anabaena_cylindrospermum = mean(anabaena_cylindrospermum),
                   other_nfixers = mean(other_nfixers)) %>% #Bare is not a living species
  mutate_if(is.numeric, log) %>%
  mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) %>%
  select(-week) 
#mutate(across(1:5, round, 0)) #green_algae:microcoleus to pull out, comment out if logging


model.2 <- list("Nweeks" = nrow(green_micro), 
                "Nspecies" = ncol(green_micro),
                "N" = green_micro #green algae and microcoleus
)

#-------------------------------------------------------------------------------------------------
#TWO SPECIES and MULTI-REACH - Gather data into STAN list format
library(abind)

#Clean and transform in a 2D dataframe
temp.spreach <- weekdata %>% 
  select(-c(1:3, 5, 8:10)) %>% 
  mutate(across(green_algae:microcoleus, round, 0)) 
# mutate(across(.cols = c("green_algae":3), .fns = log)) %>%  #logtransform
# mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) #reset the -99s

#Split data into an array by reach, then drop the reach column
spreach.array <- abind(split(temp.spreach[, -1], temp.spreach$reach), along = 3)

#Convert array into a list
spreach = plyr::alply(spreach.array,3, .dims = TRUE)


model.3 <- list("Nweeks" = nrow(spreach[["1"]]), 
                "Nreach" = length(spreach),
                "Nspecies" = ncol(spreach[["1"]]),
                "N" = spreach
)

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
                "cond" = stand_nut$cond_uS_cm
)
#-------------------------------------------------------------------------------------------------
#Run model

setwd(here::here("data_cleaning"))

options(mc.cores = parallel::detectCores())
#One year, one species, 3 reaches
# fit.m1 <-  stan(file = c("HAB_missing_weeks.stan"), data = model.1, chains = 3, iter = 10000,
#                 warmup = 5000, refresh=10, control = list(adapt_delta = 0.999,
#                                                           stepsize = 0.001,
#                                                           max_treedepth = 20))

#One year, two species, averaged reach
# fit.m2 <-  stan(file = "HAB_two_species.stan", data = model.2, chains = 3, iter = 10000,
#                 warmup = 5000, refresh=100, control = list(adapt_delta = 0.999,
#                                                            stepsize = 0.001,
#                                                            max_treedepth = 20))
# 
# #One year, two species, 3 reaches
# fit.m3 <-  stan(file = c("HAB_spreach.stan"), data = model.3, chains = 3, iter = 10000,
#                 warmup = 5000, refresh=100, control = list(adapt_delta = 0.999,
#                                                            stepsize = 0.001,
#                                                            max_treedepth = 20))

#All years, all species, averaged reach
fit.m4 <-  stan(file = "HAB_all_years.stan", data = model.4, chains = 3, iter = 10000,
                warmup = 5000, refresh=100, control = list(adapt_delta = 0.999,
                                                           stepsize = 0.001,
                                                           max_treedepth = 20))

#-------------------------------------------------------------------------------------------------
#Model checks and evaluation
library(shinystan)
library(bayesplot)
library(ggplot2)
shinystan::launch_shinystan(fit.m4)
print(fit.m4, par = "Beta")

#Save output for cleaning and visualizing in data_analysis/model_vs_real_data.R
avg.reach.fit <- fit.m2
avg.reach.output <- as.data.frame(rstan::extract(fit.m4, permuted=FALSE))


#Save raw parameter estimates
saveRDS(
  avg.reach.fit, 
  file = here::here("data/Bayes_avg_reach_fit.rds")
) 

saveRDS(
  avg.reach.output, 
  file = here::here("data/Bayes_avg_reach_output.rds")
) 
#To read: object <- readRDS(here:here("data/Bayes_avg_reach_output.rds"))
