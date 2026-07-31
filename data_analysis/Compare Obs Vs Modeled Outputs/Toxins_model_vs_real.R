###########################
#Toxins: Model vs Real Data

#This file is for creating figures that compare the collected anatoxin data with 
#modeled outputs that fill in the missing weeks per first and last year
###########################

#Package library
library(tidyverse)
library(here)
library(bayesplot)
library(ggplot2)
library(lubridate)
library(patchwork)
library(rstan)

#Read in functions
here::here("data/Functions.R")

#Read in cleaned real data 
toxinsTM <- readRDS(here::here("data/Outputs for Obs vs Real/obs_toxins_TM.rds"))
uniqueID_toxinsTM <- readRDS(here::here("data/Outputs for Sims and Model Fits/obs_toxins_data_TM.rds"))
toxinsTAC <- readRDS(here::here("data/Outputs for Obs vs Real/obs_toxins_TAC.rds"))
uniqueID_toxinsTAC <- readRDS(here::here("data/Outputs for Sims and Model Fits/obs_toxins_data_TAC.rds"))

#Read in model data (from Toxin_Models)
fit.atx.riverTM <- readRDS(here::here("data/Outputs for Obs vs Real/Anatoxin_TM_Riverwide.rds"))
fit.atx.riverTAC <- readRDS(here::here("data/Outputs for Obs vs Real/Anatoxin_TAC_Riverwide.rds"))
fit.atx.mat <- readRDS(here::here("data/Outputs for Obs vs Real/Anatoxin_Withinmat.rds"))

#Toxins from Microcoleus mats--------------------------------------------------------------
#OBSERVED DATA
obs_data_toxinsTM <- toxinsTM %>% 
  pivot_longer(cols = c("ATX_all_ug_g":"dhATXa_ug_g"), names_to = "Congener", values_to = "Concentration") %>% 
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


                                  #MODELED DATA - River-Wide

tox_params_riverTM <- as.data.frame(fit.atx.riverTM) %>% 
  dplyr::select(matches("tox_raw")) %>% 
  dplyr::mutate(across(`chain:1.tox_raw[1]`:`chain:3.tox_raw[41]`, ~ . / 1000)) %>% #backtransform for poisson
  t 

#Set up dataframe to extract week/year info from
yearweek_atxTM <- uniqueID_toxinsTM %>% 
  dplyr::rename('Total Anatoxins' = ATX_all_ug_g) %>% #Anatoxin-a = ATXa_ug_g,
                                                      #Homoanatoxin-a = HTXa_ug_g,
                                                      #Dihydroanatoxin-a = dhATXa_ug_g
  pivot_longer(cols = 3,
               names_to = "Congener", values_to = "mean") %>% 
  dplyr::mutate(time = rep(seq(41), each = length(unique(Congener))))   #41 is mat timeseries length

#Manually calculate median posteriors for microscopy proportions
tox_params2_riverTM <- as.data.frame(tox_params_riverTM) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE),
                   se_median = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  dplyr::mutate(Congener = "Total Anatoxins") %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,1])) %>%  #Change [,1] to [,2] if multiple congeners again
  left_join(yearweek_atxTM[,c("uniqueID", "Congener", "time")], by = c("Congener", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_toxinsTM[,c("year", "week", "Congener", "real_week")], 
            by = c("year", "week", "Congener")) %>% 
  arrange(time) %>% 
  mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>%
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))

#Save cleaned and summarized River-Wide toxin model latent states for visualizing 
  #predictive toxin simulations
saveRDS(tox_params2_riverTM, 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_TM_Tox_LatentStates.rds"))

#Toxins from Anabaena mats--------------------------------------------------------------
#OBSERVED DATA
obs_data_toxinsTAC <- toxinsTAC %>% 
  pivot_longer(cols = c("ATX_all_ug_g":"dhATXa_ug_g"), names_to = "Congener", values_to = "Concentration") %>% 
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


#MODELED DATA - River-Wide

tox_params_riverTAC <- as.data.frame(fit.atx.riverTAC) %>% 
  dplyr::select(matches("tox_raw")) %>% 
  dplyr::mutate(across(`chain:1.tox_raw[1]`:`chain:3.tox_raw[30]`, ~ . / 1000)) %>% #backtransform for poisson
  t 

#Set up dataframe to extract week/year info from
yearweek_atxTAC <- uniqueID_toxinsTAC %>% 
  dplyr::rename('Total Anatoxins' = ATX_all_ug_g) %>% #Anatoxin-a = ATXa_ug_g,
  #Homoanatoxin-a = HTXa_ug_g,
  #Dihydroanatoxin-a = dhATXa_ug_g
  pivot_longer(cols = 3,
               names_to = "Congener", values_to = "mean") %>% 
  dplyr::mutate(time = rep(seq(30), each = length(unique(Congener))))   #30 is tox timeseries length

#Manually calculate median posteriors for microscopy proportions
tox_params2_riverTAC <- as.data.frame(tox_params_riverTAC) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE),
                   se_median = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  dplyr::mutate(Congener = "Total Anatoxins") %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,1])) %>%  #Change [,1] to [,2] if multiple congeners again
  left_join(yearweek_atxTAC[,c("uniqueID", "Congener", "time")], by = c("Congener", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_toxinsTAC[,c("year", "week", "Congener", "real_week")], 
            by = c("year", "week", "Congener")) %>% 
  arrange(time) %>% 
  mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>%
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))

#Save cleaned and summarized River-Wide toxin model latent states for visualizing 
#predictive toxin simulations
saveRDS(tox_params2_riverTAC, 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/Riverwide_TAC_Tox_LatentStates.rds"))

                                  #MODELED DATA - Within-Mat
tox_params_mat <- as.data.frame(fit.atx.mat) %>% 
  dplyr::select(matches("tox_raw")) %>% 
  dplyr::mutate(across(`chain:1.tox_raw[1]`:`chain:3.tox_raw[41]`, ~ . / 1000)) %>% #backtransform for poisson
  t 

#Manually calculate median posteriors for microscopy proportions
tox_params2_mat <- as.data.frame(tox_params_mat) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE),
                   se_median = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  dplyr::mutate(Congener = "Total Anatoxins") %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,1])) %>%  #Change [,1] to [,2] if multiple congeners again
  left_join(yearweek_atxTM[,c("uniqueID", "Congener", "time")], by = c("Congener", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_toxinsTM[,c("year", "week", "Congener", "real_week")], 
            by = c("year", "week", "Congener")) %>% 
  arrange(time) %>% 
  mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>%
  mutate(real_week = ifelse(year == 2024, time, real_week)) %>% #manually fill in multiple skipped weeks, luckily real week = timestep in this year
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))

saveRDS(tox_params2_mat, 
        file = here::here("data/Outputs for Sims and Model Fits/Latent States/WithinMat_Tox_LatentStates.rds"))





#FIGURES--------------------------------------------------------------------------------
  
#River-Wide TM
ggplot(tox_params2_riverTM, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Congener), 
              alpha = 0.2) +
  # Latent points/lines
  geom_point(aes(colour = Congener), size = 3) +
  geom_line(aes(colour = Congener), size = 2, alpha = 0.7) +
  # Observed points/lines
  geom_point(data = subset(obs_data_toxinsTM, Congener %in% c("Total Anatoxins")), 
             aes(x = model_date, y = obs_mean, shape = Congener),
             size = 2.5) +
  geom_line(data = subset(obs_data_toxinsTM, Congener %in% c("Total Anatoxins")), 
            aes(x = model_date, y = obs_mean, group = Congener),
            size = 0.5) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(y = c(0, 40)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "Observed vs. Latent Micro Mat Tox Concentrations: River-Wide") +
  scale_colour_manual(values = c("Total Anatoxins" = "#791c55")) +
  scale_fill_manual(values = c("Total Anatoxins" = "#791c55")) +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  theme_bw()

#River-Wide TAC
ggplot(tox_params2_riverTAC, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Congener), 
              alpha = 0.2) +
  # Latent points/lines
  geom_point(aes(colour = Congener), size = 3) +
  geom_line(aes(colour = Congener), size = 2, alpha = 0.7) +
  # Observed points/lines
  geom_point(data = subset(obs_data_toxinsTAC, Congener %in% c("Total Anatoxins")), 
             aes(x = model_date, y = obs_mean, shape = Congener),
             size = 2.5) +
  geom_line(data = subset(obs_data_toxinsTAC, Congener %in% c("Total Anatoxins")), 
            aes(x = model_date, y = obs_mean, group = Congener),
            size = 0.5) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(y = c(0, 50)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "Observed vs. Latent Ana Mat Tox Concentrations: River-Wide") +
  scale_colour_manual(values = c("Total Anatoxins" = "#791c55")) +
  scale_fill_manual(values = c("Total Anatoxins" = "#791c55")) +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  theme_bw()

#Within-Mat
ggplot(tox_params2_mat, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Congener), 
              alpha = 0.2) +
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
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(y = c(0, 40)) +
  labs(x = "Date", y = "Anatoxin Concentration (ug/g)", title = "Observed vs. Latent Concentrations: Within-Mat") +
  scale_colour_manual(values = c("Total Anatoxins" = "#791c55")) +
  scale_fill_manual(values = c("Total Anatoxins" = "#791c55")) +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  theme_bw()


