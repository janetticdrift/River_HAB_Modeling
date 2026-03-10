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

#Read in environmental and microscopy data
source(here::here("data_cleaning/cleaning_HAB.R"))
#Read in anatoxin data
toxindf <- read.csv(here::here("data/HABS_anatoxins.csv"))
# #Read in gene copy data
# genesdf <- read.csv(here::here("data/qPCR_genecopies.csv"))

#Clean dataframes to feed into the toxin model

#Isolate Microcoleus mats
toxins <- toxindf %>% 
  dplyr::filter(sample_type == "Microcoleus") %>% 
  dplyr::mutate(field_date = ifelse(field_date == "7/11/23", "7/10/23", 
                                    field_date)) %>%  #Replace 7/11/23 with 7/10/23 so they are on the same week
  dplyr::filter(field_date != "6/19/24") %>% #Remove 6/19/24 as it wasn't sampled in microscopy data
  dplyr::mutate(field_date = as.Date(field_date, format="%m/%d/%y"))


anaCsplit <- toxins %>% 
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

#year 2023
year2_index <- year2 %>% 
  dplyr::mutate(timestep = dense_rank(field_date)) %>% 
  dplyr::mutate(week = timestep)

############################Note in 2023 they did not anatoxin sample the first week of %covers

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
  dplyr::select(-c(field_date, site, sample_type, timestep)) %>% 
  group_by(year) %>% 
  complete(nesting(reach), week = seq(min(week), max(week), 1L)) %>% 
  ungroup() %>%
  dplyr::mutate(reach = as.numeric(factor(reach))) 


#---------------------------------------------------------------------------------------
#CREATE MODELS

#Gather data into stan list format
anatoxin_data <- atx %>% 
  group_by(year, week) %>% 
  dplyr::summarise(ATX_all_ug_g = mean(ATX_all_ug_g, na.rm = TRUE)) %>% #Average across reaches, removing reaches where no ATX was collected
  dplyr::mutate(across(everything(), ~replace(., . == 0, 1))) %>%  #Cannot have zeros for log transforming
  replace(is.na(.), -99) %>%
  mutate(firstday = if_else(week == 1 & (year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T)


model.atx <- list("uniqueID" = nrow(anatoxin_data),
                "firstdays" = anatoxin_data$firstday,
                "Toxins" = anatoxin_data$ATX_all_ug_g,
                "Nspecies" = as.integer(ncol(matalltaxaM)-2),
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
#Run models

setwd(here::here("data_cleaning")) #Set working directory to current folder

options(mc.cores = parallel::detectCores())

#Estimate anatoxins in TM mats
fit.atx <-  stan(file = "HAB_toxins.stan", data = model.atx, chains = 3, iter = 6000,
                 warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                            max_treedepth = 15))

#Model checks and evaluation
library(shinystan)
library(bayesplot)
library(ggplot2)
library(rstantools)

#Can check posterior graphs in shinystan
shinystan::launch_shinystan(fit.atx)

