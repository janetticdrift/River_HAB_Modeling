#############
##Visualizing Modeled vs Observed Outputs##

#This file is for creating figures that compare the collected data with modeled 
#outputs that fill in the missing weeks per first and last year
#############

#Package library
library(tidyverse)
library(here)
library(bayesplot)
library(ggplot2)
library(lubridate)
library(patchwork)

#Read in functions
here::here("data/Functions.R")

#Read in observed percent cover data (obs_river_wide_viz)
source(here::here("data_cleaning/cleaning_HAB.R"))
#Read in join-matching data
alltaxatime <- readRDS(here::here("data/Outputs for Sims and Model Fits/obs_river_data.rds"))

#Read in model data (from Latent_States_Models.R)
fit.m5 <- readRDS(here::here("data/Outputs for Obs vs Real/Riverwide_Nitrate.rds"))
fit.m6 <- readRDS(here::here("data/Outputs for Obs vs Real/Riverwide_Phos.rds"))
fit.m7 <- readRDS(here::here("data/Outputs for Obs vs Real/Riverwide_Ammonium.rds"))
fit.m8 <- readRDS(here::here("data/Outputs for Obs vs Real/Riverwide_DIN.rds"))

#Model including nitrates--------------------------------------------------------------
#OBSERVED DATA - calculate average values
obs_data_all <- obs_river_wide_viz %>% 
  group_by(field_date, year, Species) %>% 
  dplyr::summarise(obs_mean = mean(Abundance), obs_SE = calcSE(Abundance)) %>% 
  ungroup() %>% 
  dplyr::group_by(year) %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - first(real_week) + 1,
         model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7)) %>% 
  dplyr::filter(Species != "bare_biofilm") %>% 
  dplyr::mutate(Species = case_when(Species == "green_algae" ~ "Green Algae",
                                    Species == "microcoleus" ~ "Microcoleus",
                                    Species == "anabaena_cylindrospermum" ~ "Anabaena",
                                    Species == "other_nfixers" ~ "Other N Fixers")) %>% 
  dplyr::mutate(Species = as.factor(Species))

#Manually calculate median posteriors for species $ cover, as well as confidence interval

#MODEL M.4 - ALL VARIABLES
params1_all <- as.data.frame(fit.m5) %>%
  dplyr::select(matches("n\\[")) %>%
  dplyr::mutate(across(1:`chain:3.n[4,45]`, exp)) %>%
  t

#Set up dataframe to extract week/year info from
yearweek <- alltaxatime %>% 
   pivot_longer(cols = c(green_algae, microcoleus, anabaena_cylindrospermum,
                         bare_biofilm, other_nfixers),
                names_to = "Species", values_to = "median") %>% 
  dplyr::mutate(time = rep(seq(45), each = length(unique(Species)))) %>% 
  dplyr::mutate(Species = case_when(Species == "green_algae" ~ "Green Algae",
                                    Species == "microcoleus" ~ "Microcoleus",
                                    Species == "anabaena_cylindrospermum" ~ "Anabaena",
                                    Species == "other_nfixers" ~ "Other N Fixers"))
 
  
params2_all_nitrate <- as.data.frame(params1_all) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE),
                   se_median = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  dplyr::mutate(Species = case_when(grepl("1,", group) ~ "Green Algae",
                             grepl("2,", group) ~ "Microcoleus",
                             grepl("3,", group) ~ "Anabaena",
                             grepl("4,", group) ~ "Other N Fixers",
                             grepl("b", group) ~ 'bare_biofilm')) %>% 
  dplyr::mutate(time = as.numeric(ifelse(grepl("b", group), str_extract(group, "[0-9]+"),
         str_extract_all(group, "[0-9]+", simplify = T)[,2]))) %>% 
  left_join(yearweek[,c("uniqueID", "Species", "time")], by = c("Species", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  dplyr::mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_all[,c("year", "week", "Species", "real_week")], by = c("year", "week", "Species")) %>% 
  arrange(time) %>% 
  dplyr::mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                              (real_week - 1) * 7 - 1, "week", week_start = 7)) %>% 
  dplyr::filter(Species != "bare_biofilm") %>% 
  dplyr::mutate(Species = as.factor(Species))


 #FIGURES--------------------------------------------------------------------------------

#Create a color palette
mycols <- c("brown", "darkolivegreen4", "darkcyan", "darkorange")
mypal <- palette(mycols)
mypal <- palette(mycols)
names(mypal) <- c("Anabaena", "Green Algae", "Microcoleus", 
                 "Other N Fixers")
colScale <- scale_color_manual(values = mypal)
filScale <- scale_fill_manual(values = mypal)

myshap <- c(16, 17, 15, 3)
names(myshap) <- c("Anabaena", "Green Algae", "Microcoleus", 
                  "Other N Fixers")
shapScale <- scale_shape_manual(values = myshap)

# Plot 1: Only show Green Algae + Other N Fixers
p1 <- ggplot(params2_all_nitrate, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(params2_all_nitrate,
                               median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Latent points/lines
  geom_point(aes(colour = Species), size = 3,
             data = transform(params2_all_nitrate, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  geom_line(aes(colour = Species), size = 2, alpha = 0.7,
            data = transform(params2_all_nitrate, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  # Observed points/lines
  geom_point(aes(y = obs_mean, shape = Species), size = 2.5,
             data = transform(obs_data_all,
                              obs_mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), obs_mean, NA))) +
  geom_line(aes(y = obs_mean, group = Species), size = 0.5,
            data = transform(obs_data_all,
                             obs_mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), obs_mean, NA))) +
  scale_y_continuous(breaks = seq(0, 200, 10)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Observed vs. Latent Abundances") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  coord_cartesian(y = c(0, 65)) +
  colScale + filScale + shapScale + theme_bw()

# Plot 2: Only show Anabaena + Microcoleus
p2 <- ggplot(params2_all_nitrate, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(params2_all_nitrate,
                               median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  geom_point(aes(colour = Species), size = 3,
             data = transform(params2_all_nitrate, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  geom_line(aes(colour = Species), size = 2, alpha = 0.7,
            data = transform(params2_all_nitrate, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  geom_point(aes(y = obs_mean, shape = Species), size = 2.5,
             data = transform(obs_data_all,
                              obs_mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), obs_mean, NA))) +
  geom_line(aes(y = obs_mean, group = Species), size = 0.5,
            data = transform(obs_data_all,
                             obs_mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), obs_mean, NA))) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  labs(x = "Date", y = "Percent Cover (%)") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  coord_cartesian(y = c(0, 25)) +
  colScale + filScale + shapScale + theme_bw()

# Combine plots and collect legends
(p1 / p2) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical") 
  
#Pulling out basic numbers: the median latent percent cover per species, per year
aggregate(median ~ Species + year, data = params2_all, max)
aggregate(obs_mean ~ Species + year, data = obs_data_all, max)

#Columns plot - need to add barebio back into params2_all for this
ggplot(obs_data_all, aes(x = model_date, y = obs_mean, fill = Species)) +
  facet_wrap(~year, scales = "free") +
  geom_col(position = "fill", width = 5) #+
  #scale_x_continuous(breaks=c(seq(1,17,2))) This was when x = week

#-------------------------------------------------------------------
  
#Latent states for phos model
#Manually calculate median posteriors for species $ cover, as well as confidence interval
params1_all <- as.data.frame(fit.m6) %>%
  dplyr::select(matches("n\\[")) %>%
  dplyr::mutate(across(1:`chain:3.n[4,45]`, exp)) %>%
  t

#Set up dataframe to extract week/year info from
yearweek <- alltaxatime %>% 
  pivot_longer(cols = c(green_algae, microcoleus, anabaena_cylindrospermum,
                        bare_biofilm, other_nfixers),
               names_to = "Species", values_to = "median") %>% 
  dplyr::mutate(time = rep(seq(45), each = length(unique(Species)))) %>% 
  dplyr::mutate(Species = case_when(Species == "green_algae" ~ "Green Algae",
                                    Species == "microcoleus" ~ "Microcoleus",
                                    Species == "anabaena_cylindrospermum" ~ "Anabaena",
                                    Species == "other_nfixers" ~ "Other N Fixers"))


params2_all_phos <- as.data.frame(params1_all) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE),
                   se_median = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  dplyr::mutate(Species = case_when(grepl("1,", group) ~ "Green Algae",
                                    grepl("2,", group) ~ "Microcoleus",
                                    grepl("3,", group) ~ "Anabaena",
                                    grepl("4,", group) ~ "Other N Fixers",
                                    grepl("b", group) ~ 'bare_biofilm')) %>% 
  dplyr::mutate(time = as.numeric(ifelse(grepl("b", group), str_extract(group, "[0-9]+"),
                                         str_extract_all(group, "[0-9]+", simplify = T)[,2]))) %>% 
  left_join(yearweek[,c("uniqueID", "Species", "time")], by = c("Species", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  dplyr::mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_all[,c("year", "week", "Species", "real_week")], by = c("year", "week", "Species")) %>% 
  arrange(time) %>% 
  dplyr::mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7)) %>% 
  dplyr::filter(Species != "bare_biofilm") %>% 
  dplyr::mutate(Species = as.factor(Species))

# Plot 1: Only show Green Algae + Other N Fixers
p3 <- ggplot(params2_all, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  #Latent points and lines
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(params2_all,
                               median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  geom_point(aes(colour = Species), size = 3,
             data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  geom_line(aes(colour = Species), size = 2, alpha = 0.7,
            data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  # Observed points and lines
  geom_point(aes(y = obs_mean, shape = Species), size = 2.5,
             data = transform(obs_data_all,
                              obs_mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), obs_mean, NA))) +
  geom_line(aes(y = obs_mean, group = Species), size = 0.5,
            data = transform(obs_data_all,
                             obs_mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), obs_mean, NA))) +
  scale_y_continuous(breaks = seq(0, 150, 10)) +
  coord_cartesian(y = c(0, 60)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Observed vs. Latent Abundances: Only Biotic Interactions") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  colScale + filScale + shapScale + theme_bw()

# Plot 2: Only show Anabaena + Microcoleus
p4 <- ggplot(params2_all, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  #Latent points and lines
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(params2_all,
                               median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  geom_point(aes(colour = Species), size = 3,
             data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  geom_line(aes(colour = Species), size = 2, alpha = 0.7,
            data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  #Observed points and lines
  geom_point(aes(y = obs_mean, shape = Species), size = 2.5,
             data = transform(obs_data_all,
                              obs_mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), obs_mean, NA))) +
  geom_line(aes(y = obs_mean, group = Species), size = 0.5,
            data = transform(obs_data_all,
                             obs_mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), obs_mean, NA))) +
  scale_y_continuous(breaks = seq(0, 150, 10)) +
  coord_cartesian(y = c(0, 25)) +
  labs(x = "Date", y = "Percent Cover (%)") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  colScale + filScale + shapScale + theme_bw()

# Combine plots and collect legends
(p3 / p4) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical") 


#-------------------------------------------------------------------
  
#Latent states for ammonium model
#Manually calculate median posteriors for species $ cover, as well as confidence interval
params1_all <- as.data.frame(fit.m7) %>%
  dplyr::select(matches("n\\[")) %>%
  dplyr::mutate(across(1:`chain:3.n[4,45]`, exp)) %>%
  t

#Set up dataframe to extract week/year info from
yearweek <- alltaxatime %>% 
  pivot_longer(cols = c(green_algae, microcoleus, anabaena_cylindrospermum,
                        bare_biofilm, other_nfixers),
               names_to = "Species", values_to = "median") %>% 
  dplyr::mutate(time = rep(seq(45), each = length(unique(Species)))) %>% 
  dplyr::mutate(Species = case_when(Species == "green_algae" ~ "Green Algae",
                                    Species == "microcoleus" ~ "Microcoleus",
                                    Species == "anabaena_cylindrospermum" ~ "Anabaena",
                                    Species == "other_nfixers" ~ "Other N Fixers"))


params2_all_ammonium <- as.data.frame(params1_all) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE),
                   se_median = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  dplyr::mutate(Species = case_when(grepl("1,", group) ~ "Green Algae",
                                    grepl("2,", group) ~ "Microcoleus",
                                    grepl("3,", group) ~ "Anabaena",
                                    grepl("4,", group) ~ "Other N Fixers",
                                    grepl("b", group) ~ 'bare_biofilm')) %>% 
  dplyr::mutate(time = as.numeric(ifelse(grepl("b", group), str_extract(group, "[0-9]+"),
                                         str_extract_all(group, "[0-9]+", simplify = T)[,2]))) %>% 
  left_join(yearweek[,c("uniqueID", "Species", "time")], by = c("Species", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  dplyr::mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_all[,c("year", "week", "Species", "real_week")], by = c("year", "week", "Species")) %>% 
  arrange(time) %>% 
  dplyr::mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7)) %>% 
  dplyr::filter(Species != "bare_biofilm") %>% 
  dplyr::mutate(Species = as.factor(Species))

# Plot 1: Only show Green Algae + Other N Fixers
p5 <- ggplot(params2_all, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(params2_all,
                               median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Latent points/lines
  geom_point(aes(colour = Species), size = 3,
             data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  geom_line(aes(colour = Species), size = 2, alpha = 0.7,
            data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  # Observed points/lines
  geom_point(aes(y = obs_mean, shape = Species), size = 2.5,
             data = transform(obs_data_all,
                              obs_mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), obs_mean, NA))) +
  geom_line(aes(y = obs_mean, group = Species), size = 0.5,
            data = transform(obs_data_all,
                             obs_mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), obs_mean, NA))) +
  scale_y_continuous(breaks = seq(0, 150, 10)) +
  coord_cartesian(y = c(0, 70)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Observed vs. Latent Abundances: Only Abiotic Interactions") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  colScale + filScale + shapScale + theme_bw()

# Plot 2: Only show Anabaena + Microcoleus
p6 <- ggplot(params2_all, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(params2_all,
                               median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  geom_point(aes(colour = Species), size = 3,
             data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  geom_line(aes(colour = Species), size = 2, alpha = 0.7,
            data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  geom_point(aes(y = obs_mean, shape = Species), size = 2.5,
             data = transform(obs_data_all,
                              obs_mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), obs_mean, NA))) +
  geom_line(aes(y = obs_mean, group = Species), size = 0.5,
            data = transform(obs_data_all,
                             obs_mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), obs_mean, NA))) +
  scale_y_continuous(breaks = seq(0, 150, 10)) +
  coord_cartesian(y = c(0, 25)) +
  labs(x = "Date", y = "Percent Cover (%)") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  colScale + filScale + shapScale + theme_bw()

# Combine plots and collect legends
(p5 / p6) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical") 

#-------------------------------------------------------------------

#Latent states for abiotic-only model minus nutrients
#Manually calculate median posteriors for species $ cover, as well as confidence interval
params1_all <- as.data.frame(fit.m8) %>%
  dplyr::select(matches("n\\[")) %>%
  dplyr::mutate(across(1:`chain:3.n[4,45]`, exp)) %>%
  t

#Set up dataframe to extract week/year info from
yearweek <- alltaxatime %>% 
  pivot_longer(cols = c(green_algae, microcoleus, anabaena_cylindrospermum,
                        bare_biofilm, other_nfixers),
               names_to = "Species", values_to = "median") %>% 
  dplyr::mutate(time = rep(seq(45), each = length(unique(Species)))) %>% 
  dplyr::mutate(Species = case_when(Species == "green_algae" ~ "Green Algae",
                                    Species == "microcoleus" ~ "Microcoleus",
                                    Species == "anabaena_cylindrospermum" ~ "Anabaena",
                                    Species == "other_nfixers" ~ "Other N Fixers"))


params2_all_DIN <- as.data.frame(params1_all) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE),
                   se_median = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  dplyr::mutate(Species = case_when(grepl("1,", group) ~ "Green Algae",
                                    grepl("2,", group) ~ "Microcoleus",
                                    grepl("3,", group) ~ "Anabaena",
                                    grepl("4,", group) ~ "Other N Fixers",
                                    grepl("b", group) ~ 'bare_biofilm')) %>% 
  dplyr::mutate(time = as.numeric(ifelse(grepl("b", group), str_extract(group, "[0-9]+"),
                                         str_extract_all(group, "[0-9]+", simplify = T)[,2]))) %>% 
  left_join(yearweek[,c("uniqueID", "Species", "time")], by = c("Species", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  dplyr::mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_all[,c("year", "week", "Species", "real_week")], by = c("year", "week", "Species")) %>% 
  arrange(time) %>% 
  dplyr::mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>% 
  dplyr::mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7)) %>% 
  dplyr::filter(Species != "bare_biofilm") %>% 
  dplyr::mutate(Species = as.factor(Species))

# Plot 1: Only show Green Algae + Other N Fixers
p7 <- ggplot(params2_all, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(params2_all,
                               median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA),
                               CIlower = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Green Algae", "Other N Fixers"), CIupper, NA))) +
  # Latent points/lines
  geom_point(aes(colour = Species), size = 3,
             data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  geom_line(aes(colour = Species), size = 2, alpha = 0.7,
            data = transform(params2_all, median = ifelse(Species %in% c("Green Algae", "Other N Fixers"), median, NA))) +
  # Observed points/lines
  geom_point(aes(y = obs_mean, shape = Species), size = 2.5,
             data = transform(obs_data_all,
                              obs_mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), obs_mean, NA))) +
  geom_line(aes(y = obs_mean, group = Species), size = 0.5,
            data = transform(obs_data_all,
                             obs_mean = ifelse(Species %in% c("Green Algae", "Other N Fixers"), obs_mean, NA))) +
  scale_y_continuous(breaks = seq(0, 150, 10)) +
  coord_cartesian(y = c(0, 65)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Observed vs. Latent Abundances: Only Abiotic Interactions - No Nutrients") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  colScale + filScale + shapScale + theme_bw()

# Plot 2: Only show Anabaena + Microcoleus
p8 <- ggplot(params2_all, aes(x = model_date, y = median)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species),
              alpha = 0.3,
              data = transform(params2_all,
                               median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA),
                               CIlower = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIlower, NA),
                               CIupper = ifelse(Species %in% c("Anabaena", "Microcoleus"), CIupper, NA))) +
  geom_point(aes(colour = Species), size = 3,
             data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  geom_line(aes(colour = Species), size = 2, alpha = 0.7,
            data = transform(params2_all, median = ifelse(Species %in% c("Anabaena", "Microcoleus"), median, NA))) +
  geom_point(aes(y = obs_mean, shape = Species), size = 2.5,
             data = transform(obs_data_all,
                              obs_mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), obs_mean, NA))) +
  geom_line(aes(y = obs_mean, group = Species), size = 0.5,
            data = transform(obs_data_all,
                             obs_mean = ifelse(Species %in% c("Anabaena", "Microcoleus"), obs_mean, NA))) +
  scale_y_continuous(breaks = seq(0, 150, 10)) +
  coord_cartesian(y = c(0, 25)) +
  labs(x = "Date", y = "Percent Cover (%)") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  colScale + filScale + shapScale + theme_bw()

# Combine plots and collect legends
(p7 / p8) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical") 

#------------------------------------------------------------------------------

#Keeping code for old style of plot, with all species combined together


# ggplot(params2_all, aes(x = model_date, y = median)) + 
#   facet_wrap(~year, scales = "free") + 
#   geom_point(aes(colour = Species), size = 3) +
#   geom_line(aes(colour = Species), size = 2, alpha = .7) +
#   # geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`,
#   #                 fill = Species), alpha = 0.3) +
#   #geom_errorbar(aes(ymin=median-se_median, ymax=median+se_median), width=.1) + 
#   geom_point(data = obs_data_all, aes(x = model_date, y = obs_median, shape = Species),
#              size = 2.5) +
#   geom_line(data = obs_data_all, aes(x = model_date, y = obs_median, group = Species),
#             size = .5) +
#   scale_y_continuous(breaks=c(seq(0,100,10))) +
#   labs(x = "Date", y = "Percent Cover (%)", title = "Observed vs. Fitted Abundances: Only Biotic Interactions") +
#   labs(color = "Modeled", fill = "Modeled", shape = "Observed") +
#   scale_color_manual(labels = c("Anabaena", "Green Algae", "Microcoleus", 
#                                 "Other N fixers"), values = c("brown", "darkolivegreen4", 
#                                                               "darkcyan", "darkorange")) +
#   scale_fill_manual(labels = c("Anabaena", "Green Algae", "Microcoleus", 
#                                "Other N fixers"), values = c("brown", "darkolivegreen4", 
#                                                              "darkcyan", "darkorange")) +
#   scale_shape_manual(labels = c("Anabaena", "Green Algae", "Microcoleus", 
#                                 "Other N Fixers"), values = c(16, 17, 15, 3))


