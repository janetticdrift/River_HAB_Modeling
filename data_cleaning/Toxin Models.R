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

#First examine relationships between anaC genes and anatoxin quantities
pseudocount <- 1e-4 #Set small number for 0s, to not return error when log-transformed

pcr_atx <- read.csv(here::here("data/qPCR_Anatoxins.csv")) %>% 
  filter(!grepl("analysis", Notes)) %>%  #Remove samples not for analysis ("floating" and "gravel")
  rename(field_date = date) %>% 
  mutate(field_date = as.Date(field_date, format="%m/%d/%y")) %>%  #Fix date data type
  mutate(
    # # Normalized anaC gene copies (copies per ng DNA)
    # normalized_anaC = ana_C.uL / DNA_conc_ng_uL,
    # normalized_anaC_rerun = ana_C.uL_rerun / DNA_conc_rerun,
    # # Normalized nif gene copies (copies per ng DNA)
    # normalized_nif = nif.uL / DNA_conc_ng_uL,
    # normalized_nif_rerun = nif.uL_rerun / DNA_conc_rerun,
    
    # Log10 normalized values
    across(c(normalized_anaC, normalized_anaC_rerun,
             normalized_nif, normalized_nif_rerun),
           ~ log10(.x + pseudocount))
    
    # Calculate Ash-Free Dry Mass for 2024
  ) %>% 
  pivot_longer(
    cols = c(normalized_anaC, normalized_anaC_rerun,
             normalized_nif, normalized_nif_rerun),
    names_to = "gene_run",
    values_to = "value"
  ) %>%
  mutate(
    gene = case_when(
      str_detect(gene_run, "anaC") ~ "anaC",
      str_detect(gene_run, "nif")  ~ "nif"
    ),
    run_type = if_else(str_detect(gene_run, "rerun"), "rerun", "original")
  ) %>%
  select(-gene_run) %>% 
  pivot_wider(names_from = gene, values_from = value) %>% 
  mutate(ATX_all_ug_L = ATX_all_ug_g * 1000)  #Transform ATX units from ug/g to ug/L
  

#ana_C vs nif, exclude zeros
ggplot(pcr_atx, aes(x = anaC, y = nif, color = mat)) +
  facet_wrap(~year, scales = "free") +
  geom_point(aes(shape = run_type), size = 2) +
  geom_smooth(data = . %>% filter(anaC != -4 & nif != -4), # filter out 0 values
    method = "lm", se = FALSE) +
  stat_regline_equation(
    data = . %>% filter(anaC != -4 & nif != -4),
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    label.x.npc = "left",
    label.y.npc = "top") +
  labs(x = "anaC gene (copies/ng)", y = "nif gene (copies/ng)") +
  theme_bw()

#ana_C vs nif, include zeros
ggplot(pcr_atx, aes(x = anaC, y = nif, color = mat)) +
  facet_wrap(~year, scales = "free") +
  geom_point(aes(shape = run_type), size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_regline_equation(
     aes(label = paste(..rr.label..)),
     formula = y ~ x,
     label.x.npc = "middle",
     label.y.npc = "top") +
  labs(x = "anaC gene (copies/ng)", y = "nif gene (copies/ng)") +
  theme_bw()

#ana_C vs Anatoxins, include zeros
ggplot(pcr_atx, aes(x = anaC, y = ATX_all_ug_L, color = mat)) +
  facet_wrap(~year, scales = "free") +
  geom_point( size = 2) +
  geom_smooth(data = . %>% filter(anaC != -4 & ATX_all_ug_L != 0),
              method = "lm", se = FALSE) +
  stat_regline_equation(
    data = . %>% filter(anaC != -4 & ATX_all_ug_L != 0),
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    label.x.npc = "left",
    label.y.npc = "top") +
  labs(x = "anaC gene (copies/ng)", y = "Total anatoxin (ug/L)") +
  scale_colour_brewer(type = "qual") +
  theme_bw()




#Compare nif presence with epithemia
nif_ediatoms <- microscopy %>%
  rename(mat = sample_type) %>% 
  mutate(mat = case_when(
    mat == "TAC" ~ "Anabaena",
    mat == "TM" ~ "Microcoleus")) %>% 
  inner_join(pcr_atx, by = c("field_date", "year", "mat")) %>%
  select(field_date, nif, e_diatoms, year, mat)

ggplot(nif_ediatoms, aes(x = nif, y = e_diatoms, color = mat)) +
  facet_wrap(~year, scales = "free") +
  geom_point() +
  geom_smooth(data = . %>% filter(nif != -4 & e_diatoms != 0),
              method = "lm", se = FALSE) +
  stat_regline_equation(
    data = . %>% filter(nif != -4 & e_diatoms != 0),
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    label.x.npc = "left",
    label.y.npc = "top") +
  labs(x = "Log nif gene (copies/ug)", y = "Epithemia (% cover microscopy)", 
       title = "Epithemia Diatoms and nif gene") +
  scale_colour_brewer(palette = "Dark2")
  

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

mat_params2_subset <- mat_params2 %>% 
  dplyr::filter(Species %in% c("anabaena_and_cylindrospermum", "e_diatoms"))


model.anaC <- list("uniqueID" = nrow(anatoxin_data),
                "firstdays" = anatoxin_data$firstday,
                "Toxins" = anatoxin_data$ATX_all_ug_g,
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
fit.toxin <-  stan(file = "HAB_toxins.stan", data = model.anaC, chains = 3, iter = 10000,
                warmup = 3000, refresh=100, control = list(adapt_delta = 0.999,
                                                                            max_treedepth = 15))
