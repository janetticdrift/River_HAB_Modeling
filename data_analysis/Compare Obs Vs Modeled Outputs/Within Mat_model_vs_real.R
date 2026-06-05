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
library(RColorBrewer)
library(lubridate)

#Read in real microscopy data 
source(here::here("data_cleaning/cleaning_HAB.R"))
#Dataframe of interest is "microscopy"

#Read in data
obs_TM_data <- readRDS(here::here("data/Outputs for Sims and Model Fits/obs_TM_data.rds"))
matmodel_TM <- readRDS(here::here("data/Outputs for Obs vs Real/WithinMat_Micro.rds"))


#Clean dataframe of observed REAL data - This one is for TM only
obs_data_mat_TM <- microscopy %>% 
  dplyr::filter(sample_type %in% "TM") %>% 
  dplyr::arrange(field_date) %>% 
  dplyr::mutate(field_date = replace(field_date, field_date == as.Date("2022-09-06"),
                              as.Date("2022-09-08"))) %>%  #Change 2022/09/06 to 09/08 so the model can summarise correctly 
  dplyr::filter(Species %in% c("Anabaena", "Epithemia Diatoms", "Geitlerinema")) %>% #Retain only Ana, Epi, and Geit
  # This code was for relativizing microscopy abundances
  # group_by(field_date, reach) %>% 
  # dplyr::mutate(Abundance = Abundance / sum(Abundance) * 100) %>% #Divide the abundances by summed total, *100
  # dplyr::mutate(Abundance = replace_na(Abundance, 0)) %>% 
  # ungroup() %>% 
  group_by(field_date, year, Species) %>% 
  dplyr::summarise(obs_mean = mean(Abundance), obs_SE = calcSE(Abundance)) %>% 
  replace(is.na(.), 0) %>% #one date had only reach, resulting in SE = NA
  ungroup() %>% 
  dplyr::group_by(year) %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - first(real_week) + 1,
                model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))

#Clean dataframe of MODEL data
#Initial cleaning of model outputs
mat_params_TM <- as.data.frame(matmodel_TM) %>% 
  dplyr::select(matches("n\\[")) %>% 
  dplyr::mutate(across(`chain:1.n[1,1]`:`chain:3.n[3,41]`, exp)) %>%  #backtransform n
  t 
    #Make sure you anazlyze n, not n_nc: n is the reconstructed latent state, and biologically meaningful
    #n_nc is just the standardized version for model construction

#Set up dataframe to extract week/year info from
yearweekTM <- obs_TM_data %>% 
  # dplyr::rename(Anabaena = anabaena_and_cylindrospermum, 
  #               'Epithemia Diatoms' = e_diatoms,
  #               Geitlerinema = geitlerinema,
  #               'Green Algae' = green_algae, 
  #               Microcoleus = microcoleus,
  #               'Non-Epithemia Diatoms' = non_e_diatoms,
  #               Nostoc = nostoc,
  #               'Other Coccoids' = other_coccoids,
  #               Rare = rare) %>% 
  pivot_longer(cols = c(3:5),
               names_to = "Species", values_to = "mean") %>% 
  dplyr::mutate(time = rep(seq(41), each = length(unique(Species))))   #41 is mat timeseries length

#Manually calculate mean posteriors for microscopy proportions
mat_params2_TM <- as.data.frame(mat_params_TM) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(mean = mean(c_across(starts_with("V")), na.rm = TRUE),
                   se_mean = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  dplyr::mutate(Species = case_when(grepl("[1,", group, fixed=TRUE) ~ 'Anabaena',
                                    grepl("[2,", group, fixed=TRUE) ~ 'Epithemia Diatoms',
                                    grepl("[3,", group, fixed=TRUE) ~ 'Geitlerinema',
                                    grepl("[4,", group, fixed=TRUE) ~ 'Green Algae',
                                    grepl("[5,", group, fixed=TRUE) ~ 'Microcoleus',
                                    grepl("[6,", group, fixed=TRUE) ~ 'Non-Epithemia Diatoms',
                                    grepl("[7,", group, fixed=TRUE) ~ 'Nostoc',
                                    grepl("[8,", group, fixed=TRUE) ~ 'Other Coccoids',
                                    grepl("[9,", group, fixed=TRUE) ~ 'Rare')) %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,2])) %>% 
  left_join(yearweekTM[,c("uniqueID", "Species", "time")], by = c("Species", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_mat_TM[,c("year", "week", "Species", "real_week")], 
            by = c("year", "week", "Species")) %>% 
  arrange(time) %>% 
  mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>%
  mutate(real_week = ifelse(year == 2024, time, real_week)) %>% #manually fill in multiple skipped weeks, luckily real week = timestep in this year
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))

#----------------------------------------------------
#Cleaning TAC model outputs

obs_data_mat_TA <- microscopy %>% 
  dplyr::filter(sample_type %in% "TAC") %>% 
  dplyr::mutate(field_date = replace(field_date, field_date == as.Date("2022-09-06"),
                                     as.Date("2022-09-08"))) %>%  #Change 2022/09/06 to 09/08 so the model can summarise correctly 
  group_by(field_date, year, Species) %>% 
  dplyr::summarise(obs_mean = mean(Abundance), obs_SE = calcSE(Abundance)) %>% 
  replace(is.na(.), 0) %>% #one date had only reach, resulting in SE = NA
  ungroup() %>% 
  dplyr::group_by(year) %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - first(real_week) + 1,
                model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                            (real_week - 1) * 7 - 1, "week", week_start = 7))

#Clean dataframe of MODEL data
#Initial cleaning of model outputs
mat_params_TA <- as.data.frame(matmodel_TA) %>% 
  dplyr::select(matches("n\\[")) %>% 
  dplyr::mutate(across(`chain:1.n[1,1]`:`chain:3.n[9,30]`, exp)) %>%  #backtransform n
  t 

#Set up dataframe to extract week/year info from
yearweekTA <- matalltaxaA %>% 
  dplyr::rename(Anabaena = anabaena_and_cylindrospermum, 
                'Epithemia Diatoms' = e_diatoms,
                Geitlerinema = geitlerinema,
                'Green Algae' = green_algae, 
                Microcoleus = microcoleus,
                'Non-Epithemia Diatoms' = non_e_diatoms,
                Nostoc = nostoc,
                'Other Coccoids' = other_coccoids,
                Rare = rare) %>% 
  pivot_longer(cols = c(3:11),
               names_to = "Species", values_to = "mean") %>% 
  mutate(time = rep(seq(30), each = length(unique(Species))))  #30 is mat timeseries length

#Manually calculate mean posteriors for microscopy proportions
mat_params2_TA <- as.data.frame(mat_params_TA) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(mean = mean(c_across(starts_with("V")), na.rm = TRUE),
                   se_mean = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>%
  dplyr::mutate(Species = case_when(grepl("[1,", group, fixed=TRUE) ~ 'Anabaena',
                                    grepl("[2,", group, fixed=TRUE) ~ 'Epithemia Diatoms',
                                    grepl("[3,", group, fixed=TRUE) ~ 'Geitlerinema',
                                    grepl("[4,", group, fixed=TRUE) ~ 'Green Algae',
                                    grepl("[5,", group, fixed=TRUE) ~ 'Microcoleus',
                                    grepl("[6,", group, fixed=TRUE) ~ 'Non-Epithemia Diatoms',
                                    grepl("[7,", group, fixed=TRUE) ~ 'Nostoc',
                                    grepl("[8,", group, fixed=TRUE) ~ 'Other Coccoids',
                                    grepl("[9,", group, fixed=TRUE) ~ 'Rare')) %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,2])) %>% 
  left_join(yearweekTA[,c("uniqueID", "Species", "time")], by = c("Species", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_mat_TA[,c("year", "week", "Species", "real_week")], 
            by = c("year", "week", "Species")) %>% 
  arrange(time) %>% 
  mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>%
  mutate(real_week = ifelse(year == 2024, time, real_week)) %>% #manually fill in multiple skipped weeks, luckily real week = timestep in this year
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))

#FIGURES--------------------------------------------------------------------------------

#Observed Data plots
###Separate out reaches for OBSERVED data
pivot_matTM <- yearmatdata_TM %>% 
  pivot_longer(cols = c(Anabaena:Rare),
               names_to = "Species", values_to = "Abundance") %>% 
  na.omit()

ggplot(subset(pivot_matTM, Species %in% c("Microcoleus")),
              aes(x = week, y = Abundance))+
  facet_grid(reach~year) +
  geom_point(aes(colour = Species)) +
  geom_line(aes(colour = Species))

model <- aov(Microcoleus ~ reach, data = subset(yearmatdata_TM, reach %in% c("1S", "3")))
summary(model)


#Model data plots

#Create a color palette
mycols <- c("brown", "#07707B", "#F6926A")
mypal <- palette(mycols)
names(mypal) <- c("Anabaena", "Epithemia Diatoms", "Geitlerinema")
colScale <- scale_color_manual(values = mypal)
filScale <- scale_fill_manual(values = mypal)

###TM: Ana + Geit + Epithemia
ggplot(mat_params2_TM, aes(x = model_date, y = mean)) +
  facet_wrap(~year, scales = "free_x") +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, fill = Species), alpha = 0.2) +
  # Latent points/lines
  geom_point(aes(colour = Species), size = 3) +
  geom_line(aes(colour = Species), size = 2, alpha = 0.7) +
  # Observed points/lines
  geom_point(data = obs_data_mat_TM, 
             aes(x = model_date, y = obs_mean, shape = Species), #shape = Species in aes
             size = 2.5) +
  geom_line(data = obs_data_mat_TM, 
            aes(x = model_date, y = obs_mean, group = Species),
            size = 0.5) +
  scale_y_continuous(breaks = seq(0, 100, 10)) +
  coord_cartesian(ylim = c(0,20)) +
  labs(x = "Date", y = "Relative Abundance (%)", title = "Observed vs. Latent Abundances") +
  labs(color = "Latent", fill = "Latent", shape = "Observed") +
  colScale + filScale + theme_bw()


###First 6 species TAC---------------------------

ggplot(subset(mat_params2_TA, Species %in% c("anabaena_and_cylindrospermum", 
                                          "e_diatoms", "geitlerinema", "green_algae",
                                          "leptolyngbya", "microcoleus")),
       aes(x = model_date, y = mean)) + 
  facet_wrap(~year, scales = "free") + 
  geom_point(aes(colour = Species), size = 3) +
  geom_line(aes(colour = Species), size = 2, alpha = .7) +
  geom_errorbar(aes(ymin=mean-se_mean, ymax=mean+se_mean), width=.1) +
  geom_point(data = subset(obs_data_mat_TA, Species %in% c("anabaena_and_cylindrospermum", 
                                                           "e_diatoms", "geitlerinema", "green_algae",
                                                           "leptolyngbya", "microcoleus")), 
             aes(x = model_date, y = obs_mean, shape = Species), #shape = Species in aes
             size = 2.5) +
  geom_line(data = subset(obs_data_mat_TA, Species %in% c("anabaena_and_cylindrospermum", 
                                                          "e_diatoms", "geitlerinema", "green_algae",
                                                          "leptolyngbya", "microcoleus")),
            aes(x = model_date, y = obs_mean, group = Species),
            size = .5) +
  scale_y_continuous(breaks=c(seq(0,100,10))) +
  labs(x = "Date", y = "Proportion (%)", title = "Target Microcoleus:Averaged Reaches") +
  labs(color = "Modeled", fill = "Modeled", shape = "Observed") +
  scale_colour_brewer(labels = c("Anabaena", "Epithemia", "Geitlerinema", 
                                 "Green Algae", "Leptolyngbya", "Microcoleus",
                                 "Non-Epithemia", "Nostoc", "Oscillatoria", "Other Coccoids",
                                 "Rare"), palette = "Set3") +
  scale_shape_manual(labels = c("Anabaena", "Epithemia", "Geitlerinema", 
                                "Green Algae", "Leptolyngbya", "Microcoleus",
                                "Non-Epithemia", "Nostoc", "Oscillatoria", "Other Coccoids",
                                "Rare"), values = c(16, 17, 15, 3, 5, 10)) +
  coord_cartesian(ylim = c(0,16)) 

