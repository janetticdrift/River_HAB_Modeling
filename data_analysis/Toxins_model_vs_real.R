###########################
#Toxins: Model vs Real Data

#This file is for creating figures that compare the collected anatoxin data with 
#modeled outputs that fill in the missing weeks per first and last year
###########################

#Package library
library(tidyverse)
library(here)
library(shinystan)
library(bayesplot)
library(ggplot2)
library(lubridate)
library(patchwork)
library(rstan)

#Read in functions
here::here("data/Functions.R")

#Read in real data (from Toxin_Models.R - should be put in other files eventually)
toxins

# #Read in real microscopy data 
# source(here::here("data_cleaning/cleaning_HAB.R"))
# #Dataframe of interest is "toxins"

#Read in model data (from Missing Week Estimates)
fit.atx1 <- rstan::extract(fit.atx, permuted=FALSE)
#readRDS(here::here("data/Anatoxin_AllVariables.rds"))

#MODEL INCLUDING JUST ATX_ALL--------------------------------------------------------------
#OBSERVED DATA
obs_data_toxins <- toxins %>% 
  pivot_longer(cols = c(6:9), names_to = "Congener", values_to = "Concentration") %>% 
  group_by(field_date, year, sample_type, Congener) %>% 
  dplyr::summarise(obs_mean = mean(Concentration), obs_SE = calcSE(Concentration)) %>% 
  ungroup() %>% 
  dplyr::group_by(year) %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - first(real_week) + 1,
                model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7)) %>% 
  dplyr::mutate(Congener = case_when(Congener == "ATX_all_ug_g" ~ "Total Anatoxins",
                                    Congener == "ATXa_ug_g" ~ "Anatoxin-a",
                                    Congener == "HTXa_ug_g" ~ "Homoanatoxin-a",
                                    Congener == "dhATXa_ug_g" ~ "Dihydroanatoxin-a")) %>% 
  dplyr::mutate(Congener = as.factor(Congener))


#Clean dataframe of MODELED data
#Initial cleaning of model outputs
tox_params_TM <- as.data.frame(fit.atx1) %>% 
  dplyr::select(matches("tox_raw")) %>% 
  #dplyr::mutate(across(`chain:1.tox_raw[1]`:`chain:3.tox_raw[41]`, ~ . / 1000)) %>%
  #dplyr::mutate(across(`chain:1.tox_raw[1]`:`chain:3.tox_raw[41]`, exp)) %>%  #backtransform tox
  t 
#Make sure you anazlyze tox, not tox_nc: tox is the reconstructed latent state, and biologically meaningful
#tox_nc is just the standardized version for model construction

#Set up dataframe to extract week/year info from
yearweek_atx <- anatoxin_data %>% 
  dplyr::rename('Total Anatoxins' = ATX_all_ug_g) %>% #Anatoxin-a = ATXa_ug_g,
                                                      #Homoanatoxin-a = HTXa_ug_g,
                                                      #Dihydroanatoxin-a = dhATXa_ug_g
  pivot_longer(cols = 3,
               names_to = "Congener", values_to = "mean") %>% 
  dplyr::mutate(time = rep(seq(41), each = length(unique(Congener))))   #41 is mat timeseries length

#Manually calculate mean posteriors for microscopy proportions
tox_params2_TM <- as.data.frame(tox_params_TM) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(mean = mean(c_across(starts_with("V")), na.rm = TRUE),
                   se_mean = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  dplyr::mutate(Congener = "Total Anatoxins") %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,1])) %>%  #Change [,1] to [,2] if multiple congeners again
  left_join(yearweek_atx[,c("uniqueID", "Congener", "time")], by = c("Congener", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_toxins[,c("year", "week", "Congener", "real_week")], 
            by = c("year", "week", "Congener")) %>% 
  arrange(time) %>% 
  mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>%
  mutate(real_week = ifelse(year == 2024, time, real_week)) %>% #manually fill in multiple skipped weeks, luckily real week = timestep in this year
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))


#FIGURES--------------------------------------------------------------------------------

mycols <- c("brown", "darkolivegreen4", "darkorange", "chartreuse3")
mypal <- palette(mycols)
names(mypal) = c("Total Anatoxins", "Anatoxin-a", "Homoanatoxin-a", 
                 "Dihydroanatoxin-a")
colScale <- scale_color_manual(values = mypal)
filScale <- scale_fill_manual(values = mypal)

myshap <- c(16, 17, 15, 3)
names(myshap) = c("Total Anatoxins", "Anatoxin-a", "Homoanatoxin-a", 
                  "Dihydroanatoxin-a")
shapScale <- scale_shape_manual(values = myshap)

###Anatoxins
ggplot(tox_params2_TM, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Congener), alpha = 0.2) +
  # Latent points/lines
  geom_point(aes(colour = Congener), size = 3) +
  geom_line(aes(colour = Congener), size = 2, alpha = 0.7) +
  # Observed points/lines
  geom_point(data = subset(obs_data_toxins, Congener %in% c("Total Anatoxins")), 
             aes(x = model_date, y = obs_mean, shape = Congener),
             size = 2.5) +
  geom_line(data = subset(obs_data_toxins, Congener %in% c("Total Anatoxins")), 
            aes(x = model_date, y = obs_mean, group = Congener),
            size = 0.5) +
  scale_y_continuous(breaks = seq(0, 200, 10)) +
  labs(x = "Date", y = "Anatoxin Concentration", title = "Poisson Observed vs. Latent Concentrations") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  colScale + filScale + shapScale

