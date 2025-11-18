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

#Clean dataframe of observed REAL data - This one is for TM only
obs_data_mat_TM <- microscopy %>% 
  dplyr::filter(sample_type %in% "TM") %>% 
  pivot_longer(cols = c(anabaena_and_cylindrospermum:rare),
               names_to = "Species", values_to = "Abundance") %>% 
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
mat_params <- as.data.frame(rstan::extract(fit.m1, permuted=FALSE)) %>% 
  dplyr::select(-c(1:`chain:3.Beta[11,11]`)) %>%
  dplyr::select(-c(`chain:1.lp__`:`chain:3.lp__`)) %>% 
  dplyr::mutate(across(`chain:1.n[1,1]`:`chain:3.n[11,26]`, exp)) %>%  #backtransform n
  t 
    #Make sure you anazlyze n, not n_nc: n is the reconstructed latent state, and biologically meaningful
    #n_nc is just the standardized version for model construction

#Set up dataframe to extract week/year info from
yearweek <- matalltaxaM %>% 
  pivot_longer(cols = c(anabaena_and_cylindrospermum:rare),
               names_to = "Species", values_to = "mean") %>% 
  mutate(time = rep(seq(26), each = length(unique(Species))))

#Manually calculate mean posteriors for microscopy proportions
mat_params2 <- as.data.frame(mat_params) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(mean = mean(c_across(starts_with("V")), na.rm = TRUE),
                   se_mean = calcSE(c_across(starts_with("V")))) %>% 
  dplyr::mutate(Species = case_when(grepl("[1,", group, fixed=TRUE) ~ 'anabaena_and_cylindrospermum',
                             grepl("[2,", group, fixed=TRUE) ~ 'e_diatoms',
                             grepl("[3,", group, fixed=TRUE) ~ 'geitlerinema',
                             grepl("[4,", group, fixed=TRUE) ~ 'green_algae',
                             grepl("[5,", group, fixed=TRUE) ~ 'leptolyngbya',
                             grepl("[6,", group, fixed=TRUE) ~ 'microcoleus',
                             grepl("[7,", group, fixed=TRUE) ~ 'non_e_diatoms',
                             grepl("[8,", group, fixed=TRUE) ~ 'nostoc',
                             grepl("[9,", group, fixed=TRUE) ~ 'oscillatoria',
                             grepl("[10,", group, fixed=TRUE) ~ 'other_coccoids',
                             grepl("[11,", group, fixed=TRUE) ~ 'rare')) %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,2])) %>% 
  left_join(yearweek[,c("uniqueID", "Species", "time")], by = c("Species", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_mat_TM[,c("year", "week", "Species", "real_week")], 
            by = c("year", "week", "Species")) %>% 
  arrange(time) %>% 
  mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))


#FIGURES--------------------------------------------------------------------------------

###First 6 species
ggplot(subset(mat_params2, Species %in% c("anabaena_and_cylindrospermum", 
                                       "e_diatoms", "geitlerinema", "green_algae",
                                       "leptolyngbya", "microcoleus")),
              aes(x = model_date, y = mean)) + 
  facet_wrap(~year, scales = "free") + 
  geom_point(aes(colour = Species), size = 3) +
  geom_line(aes(colour = Species), size = 2, alpha = .7) +
  geom_errorbar(aes(ymin=mean-se_mean, ymax=mean+se_mean), width=.1) +
  geom_point(data = subset(obs_data_mat_TM, Species %in% c("anabaena_and_cylindrospermum", 
                                                           "e_diatoms", "geitlerinema", "green_algae",
                                                           "leptolyngbya", "microcoleus")), 
                           aes(x = model_date, y = obs_mean, shape = Species), #shape = Species in aes
             size = 2.5) +
  geom_line(data = subset(obs_data_mat_TM, Species %in% c("anabaena_and_cylindrospermum", 
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


###Last 5 species#########################################################################

ggplot(subset(mat_params2, Species %in% c("non_e_diatoms", "nostoc", "oscillatoria",
                                          "other_coccoids", "rare")),
       aes(x = model_date, y = mean)) + 
  facet_wrap(~year, scales = "free") + 
  geom_point(aes(colour = Species), size = 3) +
  geom_line(aes(colour = Species), size = 2, alpha = .7) +
  geom_errorbar(aes(ymin=mean-se_mean, ymax=mean+se_mean), width=.1) +
  geom_point(data = subset(obs_data_mat_TM, Species %in% c("non_e_diatoms", "nostoc", "oscillatoria",
                                                           "other_coccoids", "rare")), 
             aes(x = model_date, y = obs_mean, shape = Species), 
             size = 2.5) +
  geom_line(data = subset(obs_data_mat_TM, Species %in% c("non_e_diatoms", "nostoc", "oscillatoria",
                                                          "other_coccoids", "rare")),
            aes(x = model_date, y = obs_mean, group = Species),
            size = .5) +
  scale_y_continuous(breaks=c(seq(0,100,10))) +
  labs(x = "Date", y = "Proportion (%)", title = "Target Microcoleus:Averaged Reaches") +
  labs(color = "Modeled", fill = "Modeled", shape = "Observed") +
  scale_colour_brewer(labels = c("Non-Epithemia", "Nostoc", "Oscillatoria", "Other Coccoids",
                                 "Rare"), palette = "Set3") +
  scale_shape_manual(labels = c("Non-Epithemia", "Nostoc", "Oscillatoria", "Other Coccoids",
                                "Rare"), values = c(16, 17, 15, 3, 5, 10)) +
  coord_cartesian(ylim = c(0,10)) 
