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
pseudocount <- 0.00001

#Gather data into stan list format
anatoxin_data <- atx %>% 
  group_by(year, week) %>% 
  dplyr::summarise(ATX_all_ug_g = mean(ATX_all_ug_g, na.rm = TRUE)) %>% #Average across reaches, removing reaches where no ATX was collected
  dplyr::mutate(ATX_all_ug_g = round(ATX_all_ug_g, digits = 3), #Editing data for poisson
                ATX_all_ug_g = ATX_all_ug_g*1000) %>%
  # dplyr::mutate(across(ATX_all_ug_g,
  #                      ~ . + pseudocount)) %>% #Cannot have zeros for log transforming
  # dplyr::mutate(across(ATX_all_ug_g, log)) %>%
  dplyr::mutate(across(everything(), ~replace(.x, is.nan(.x), -99))) %>% 
  mutate(firstday = if_else(week == 1 & (year == 2023 | year == 2024), 1, 0)) %>% 
  relocate(firstday) %>% 
  unite("uniqueID", c(year, week), sep = "_", remove=T) %>% 
  dplyr::mutate(is_obs  = ifelse(ATX_all_ug_g == -99, 0, 1), #Editing data for poisson
                ATX_all_ug_g = ifelse(ATX_all_ug_g == -99, 0, ATX_all_ug_g))


#Gather latent states of microscopy abundances
matmodel_TM <- readRDS(here::here("data/WithinMat_Micro.rds"))

TM_latent1 <- as.data.frame(matmodel_TM) %>% 
  dplyr::select(matches("n\\[")) %>% 
  dplyr::mutate(across(`chain:1.n[1,1]`:`chain:3.n[3,41]`, exp)) %>%  #backtransform n
  t 

TM_latent <- as.data.frame(TM_latent1) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(mean = mean(c_across(starts_with("V")), na.rm = TRUE)) %>% 
  dplyr::mutate(Species = case_when(grepl("[1,", group, fixed=TRUE) ~ 'Anabaena',
                                    grepl("[2,", group, fixed=TRUE) ~ 'Epithemia Diatoms',
                                    grepl("[3,", group, fixed=TRUE) ~ 'Geitlerinema')) %>% 
  mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,2])) %>% 
  dplyr::select(-group) %>% 
  dplyr::mutate(mean = log(mean)) %>% 
  pivot_wider(names_from = Species, values_from = mean) %>% 
  arrange(time)

#Create design matrix
X <- cbind(
  intercept = 1,
  TM_latent[, -1],  #Abundances are log-transformed
  nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)],
  phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)],
  ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)],
  discharge = discharge$stand_discharge[-c(14:15, 29:30)],
  temp = stand_nut$temp_C[-c(14:15, 29:30)],
  cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)],
  rad = swradiation$stand_rad[-c(14:15, 29:30)]
)

#Combine with other information into model list
model.atx <- list("uniqueID" = nrow(anatoxin_data),
                  "is_obs" = anatoxin_data$is_obs, #poisson edit
                  "firstdays" = anatoxin_data$firstday,
                  "Toxins" = as.integer(anatoxin_data$ATX_all_ug_g), #poisson edit needs as.integer
                  "Nspecies" = as.integer(ncol(matalltaxaM)-2),
                  "X" = X,
                  "Npredictors" = ncol(X)
)

# model.atx <- list("uniqueID" = nrow(anatoxin_data),
#                  # "is_obs" = anatoxin_data$is_obs,
#                 "firstdays" = anatoxin_data$firstday,
#                 "Toxins" = anatoxin_data$ATX_all_ug_g, #must use as.integer if poisson
#                 "Nspecies" = as.integer(ncol(matalltaxaM)-2),
#                 "N" = TM_latent[,-(1)],
#                 "nitrate" = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)], #Can subset 2024 out with 29:45
#                 "phos" = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)], #and also first two weeks of 2023 and 2024
#                 "ammonium" = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)], #Which is 14:15 and 29:30
#                 "discharge" = discharge$stand_discharge[-c(14:15, 29:30)],
#                 "temp" = stand_nut$temp_C[-c(14:15, 29:30)],
#                 "cond" = stand_nut$cond_uS_cm[-c(14:15, 29:30)],
#                 "rad" = swradiation$stand_rad[-c(14:15, 29:30)]
# )


#-------------------------------------------------------------------------------------------------
#Run models

setwd(here::here("data_cleaning")) #Set working directory to current folder

options(mc.cores = parallel::detectCores())

#Set starting values
init_fun_atx <- function() list(
  sigma_p = 0.5,    
  #sigma_o = 0.5, #not needed for poisson
  Beta0 = 0,
  Beta1 = 0,     # small start for species abundances
  Beta2 = 0,
  Beta3 = 0,
  tox_nc = rep(0, nrow(anatoxin_data)) #nrow is the time length, 4.79 is the mean anatoxin concentration
 )

#Estimate anatoxins in TM mats
fit.atx <-  stan(file = "HAB_toxins_poisson.stan", data = model.atx, chains = 3, iter = 6000,
                 warmup = 3000, refresh=100, init = init_fun_atx, control = list(adapt_delta = 0.999,
                                                            max_treedepth = 15))

#Model checks and evaluation
library(shinystan)
library(bayesplot)
library(ggplot2)
library(rstantools)

#Can check posterior graphs in shinystan
shinystan::launch_shinystan(as.shinystan(fit.atx))


mcmc_intervals(
  as.array(fit.atx),
  pars = c("Ntheta", "Ptheta", "Atheta")) 

mcmc_intervals(
  as.array(fit.atx),
  pars = c("Dtheta", "Ttheta", "Ctheta", "Rtheta")) 

mcmc_intervals(
  as.array(fit.atx),
  pars = c("Beta1", "Beta2", "Beta3") )

mcmc_intervals(
  as.array(fit.atx),
  pars = c("BetaAna", "BetaEpi", "BetaGeit") )


mcmc_intervals(
  as.array(fit.atx),
  pars = c("sigma_p", "sigma_o") )



#For building the observation vs latent state plots
saveRDS(rstan::extract(fit.atx, permuted=FALSE), 
        file = here::here("data/Anatoxin_AllVariables.rds"))
#For building the latent state vs predictions plots
saveRDS(rstan::extract(fit.atx), 
        file = here::here("data/Anatoxin_AllVar_predictions.rds"))



##### Exploratory business#####
lag_df <- data.frame(time = TM_latent$time,
                     ATX = exp(anatoxin_data$ATX_all_ug_g), 
                     Anabaena_mat = exp(TM_latent$Anabaena),
                     Anabaena_river = exp(alltaxatime$anabaena_cylindrospermum[-c(14:15, 29:30)]))

#Plots
ccf(lag_df$Anabaena_mat, lag_df$ATX, lag.max = 5)
ccf(lag_df$Anabaena_river, lag_df$ATX, lag.max = 5)

#Create lags
lag_df <- lag_df %>%
  dplyr::mutate(Ana_lag1 = lag(Anabaena_river, 1),
                Ana_lag2 = lag(Anabaena_river, 2)) %>% 
  dplyr::select(!Anabaena_mat)

lag_dfpivot <- lag_df %>% 
  pivot_longer(cols = starts_with("Ana"),
  names_to = "lag",
  values_to = "Abundance"
)

#Plot regression
ggplot(lag_dfpivot, aes(x = Abundance, y = ATX, color = lag)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Abundance (lagged)", y = "Toxin") +
  theme_bw()

###
#Examine DeltaATX ~ DeltaPercent + DeltaEnv

obs_TM_mat <- microscopy %>% 
  dplyr::filter(sample_type %in% "TM") %>% 
  dplyr::arrange(field_date) %>% 
  dplyr::mutate(field_date = replace(field_date, field_date == as.Date("2022-09-06"),
                                     as.Date("2022-09-08"))) %>%  #Change 2022/09/06 to 09/08 so the model can summarise correctly 
  dplyr::filter(Species %in% c("Anabaena", "Epithemia Diatoms", "Geitlerinema")) %>% 
  pivot_wider(names_from = Species, values_from = Abundance)

common_rows <- inner_join(obs_TM_mat, toxins, by = c("field_date", "reach", "year")) %>% 
  dplyr::slice(-17) #Remove odd duplicate row

deltadf <- data.frame(date = common_rows$field_date,
                      year = common_rows$year,
                      reach = common_rows$reach,
                      ATX = common_rows$ATX_all_ug_g, 
                      Anabaena = common_rows$Anabaena,
                      Epithemia = common_rows$`Epithemia Diatoms`,
                      Geitlerinema = common_rows$Geitlerinema) %>% 
  group_by(year, reach) %>% 
  arrange(date) %>% 
  dplyr::mutate(across(ATX:Geitlerinema, ~ .x - lag(.x), .names = "change_{.col}")) %>% 
  pivot_longer(c(9:11), names_to = "Species", values_to = "changepercent")


ggplot(deltadf, aes(x = changepercent, y = change_ATX, group = Species, color = Species)) +
  facet_grid(reach~year, scales = "free_x") +
  geom_point(size = 2) +
  #geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F) +
  labs(x = "Change in Within-Mat Percent Abundance", y = "Change in ATX") +
  coord_cartesian(ylim = c(-120, 120))

common_rows <- inner_join(obs_TM_mat, toxins, by = c("field_date", "reach", "year")) %>% 
  dplyr::slice(-17) #Remove odd duplicate row

deltadf <- data.frame(date = common_rows$field_date,
                      year = common_rows$year,
                      reach = common_rows$reach,
                      ATX = common_rows$ATX_all_ug_g, 
                      Anabaena = common_rows$Anabaena,
                      Epithemia = common_rows$`Epithemia Diatoms`,
                      Geitlerinema = common_rows$Geitlerinema) %>% 
  group_by(year, reach) %>% 
  arrange(date) %>% 
  dplyr::mutate(across(ATX:Geitlerinema, ~ .x - lag(.x), .names = "change_{.col}")) %>% 
  pivot_longer(c(9:11), names_to = "Species", values_to = "changepercent")
  
                      
