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

#Read in real microscopy data 
source(here::here("data_cleaning/cleaning_HAB.R"))
#Dataframe of interest is "microscopy"

#Read in model data (from Missing Week Estimates)
#Microcoleus
fit.m1 <- readRDS(here::here("data/Bayes_all_year_fit.rds")) #averaged reaches
fit.m1.1S <- readRDS(here::here("data/Bayes_all_year_fit.rds")) #reach 1S
fit.m1.3 <- readRDS(here::here("data/Bayes_all_year_fit.rds")) #reach 3
fit.m1.4 <- readRDS(here::here("data/Bayes_all_year_fit.rds")) #reach 4
#Anabaena
fit.m2 <- readRDS(here::here("data/Bayes_all_year_fit.rds")) #averaged reaches
fit.m2.1S <- readRDS(here::here("data/Bayes_all_year_fit.rds")) #reach 1S
fit.m2.3 <- readRDS(here::here("data/Bayes_all_year_fit.rds")) #reach 3
fit.m2.4 <- readRDS(here::here("data/Bayes_all_year_fit.rds")) #reach 4

#Read in join-matching data (from Missing Week Estimates)
matalltaxaM <- readRDS(here::here("data/matalltaxaM.rds"))

#Clean dataframe of observed data
obs_data_mat <- microscopy %>% 
  group_by(field_date, year, Species) %>% 
  dplyr::summarise(obs_mean = mean(Abundance), obs_SE = calcSE(Abundance)) %>% 
  ungroup() %>% 
  dplyr::group_by(year) %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - first(real_week) + 1,
                model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7)) %>% 
  dplyr::filter(Species != "bare_biofilm")