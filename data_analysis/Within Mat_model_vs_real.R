#############
##Visualizing Modeled vs Observed Outputs--Within-Mat##

#This file is for creating figures that compare the collected microscopy data with modeled 
#outputs
#############

#Package library
library(tidyverse)
library(here)
library(shinystan)
library(bayesplot)
library(ggplot2)
library(lubridate)

#Necessary functions
#SE function
calcSE <- function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#Read in real microscopy data 
source(here::here("data_cleaning/cleaning_HAB.R"))
#Dataframe of interest is "microscopy"

#Run models (from Missing Week Estimates)
#And read in join-matching data (matalltaxaM)
source(here::here("data_cleaning/Missing Week Estimates.R"))

#Clean dataframe of observed REAL data
obs_data_mat <- microscopy %>% 
  pivot_longer(cols = c(anabaena_and_cylindrospermum:rare),
               names_to = "Species", values_to = "Abundance") %>% 
  group_by(field_date, year, Species) %>% 
  dplyr::summarise(obs_mean = mean(Abundance), obs_SE = calcSE(Abundance)) %>% 
  ungroup() %>% 
  dplyr::group_by(year) %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - first(real_week) + 1,
                model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))

#Clean dataframe of MODEL data
#Manually calculate mean posteriors for species $ cover, as well as confidence interval
params1_all <- as.data.frame(rstan::extract(fit.m4, permuted=FALSE)) %>% 
  dplyr::select(-c(1:`chain:3.Beta_off[4,4]`)) %>% 
  dplyr::select(-c(`chain:1.lp__`:`chain:3.lp__`)) %>% 
  dplyr::select(-c(`chain:1.Ntheta[1]`:`chain:3.Beta[4,4]`)) %>% 
  dplyr::mutate(across(1:`chain:3.n[4,45]`, exp)) %>% 
  t 
