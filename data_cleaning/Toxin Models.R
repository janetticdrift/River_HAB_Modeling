########################################
###Anatoxin Dataframes for STAN Model###
########################################
#This file contains code reading in and cleaning data on anatoxins, gene count, microscopy data,
#and environmental covariates to predict anatoxins

#Packages----
library(rstan)
library(tidyverse)
library(dataRetrieval)

#Read in environmental and microscopy data
source(here::here("data_cleaning/cleaning_HAB.R"))
#Read in anatoxin data
toxindf <- read.csv(here::here("data/HABS_anatoxins.csv"))
#Read in gene copy data
genesdf <- read.csv(here::here("data/qPCR_genecopies.csv"))

#Clean dataframes to feed into the toxin model

toxins <- toxindf %>% 
  dplyr::filter(sample_type == "TM") %>% 
  dplyr::mutate(field_date = as.Date(field_date)) 
  

anaCsplit <- toxins %>% 
  group_split(year)
year1 <- anaCsplit[[1]]
year2 <- anaCsplit[[2]]
year3 <- anaCsplit[[3]]

week_from_step <- function(x) {
  case_when(
    timestep == 1 ~ 1,
    timestep == 2 ~ 3,
    timestep == 3 ~ 5,
    timestep == 4 ~ 7,
    timestep == 5 ~ 9,
    timestep == 6 ~ 11,
    timestep == 7 ~ 13,
    timestep == 8 ~ 15
  )
}

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

anaC <- rbind(year1_index, year2_index, year3_index) %>% 
  dplyr::select(-c(field_date, site, sample_type)) %>% 
  group_by(year) %>% 
  complete(nesting(reach), week = seq(min(week), max(week), 1L)) %>% 
  #replace(is.na(.), -99) %>% 
  ungroup() %>%
  mutate(reach = as.numeric(factor(reach)))

#I needed to summarise before the -99s were added

#---------------------------------------------------------------------------------------
#CREATE MODELS


anatoxin_data <- anaC %>% 
  group_by(year, week) %>% 
  dplyr::summarise(ATX_all_ug_g = mean(ATX_all_ug_g, na.rm = TRUE)) %>% #Average across reaches
  replace(is.na(.), -99) %>%
  mutate(firstday = if_else(week == 1 & (year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T)


model.anaC <- list("uniqueID" = nrow(anatoxin_data),
                "firstdays" = anatoxin_data$firstday,
                "Toxins" = anatoxin_data[,-(1:2)],
                "nitrate" = stand_nut$nitrate_mg_N_L,
                "phos" = stand_nut$oPhos_ug_P_L,
                "ammonium" = stand_nut$ammonium_mg_N_L,
                "discharge" = discharge$stand_discharge,
                "temp" = stand_nut$temp_C,
                "cond" = stand_nut$cond_uS_cm,
                "rad" = swradiation$stand_rad
)


#-------------------------------------------------------------------------------------------------
#Run models

setwd(here::here("data_cleaning")) #Set working directory to current folder

options(mc.cores = parallel::detectCores())

#Estimate anatoxins in TM mats
fit.m1 <-  stan(file = "HAB_toxins.stan", data = model.1, chains = 3, iter = 10000,
                warmup = 3000, refresh=100, init = init_fun, control = list(adapt_delta = 0.999,
                                                                            max_treedepth = 15))