#############
##Visualizing Modeled vs Observed Outputs##

#This file is for creating figures that compare the collected data with modeled 
#outputs that fill in the missing weeks per first and last year
#############

#Package library
library(tidyverse)
library(here)
library(shinystan)
library(bayesplot)
library(ggplot2)
library(lubridate)

#Read in real data (from allcoverdataplot dataset in SFE_raw_analysis)
coverpercent <- readRDS(here::here("data/allcoverdataplot.rds"))

#Read in model data (from Missing Week Estimates)
fit.m2 <- readRDS(here::here("data/Bayes_avg_reach_fit.rds"))
fit.m4 <- readRDS(here::here("data/Bayes_all_year_fit.rds"))


#Read in join-matching data (from Missing Week Estimates)
alltaxatime <- readRDS(here::here("data/alltaxatime.rds"))

#Necessary functions
#SE function
calcSE <- function(x){
  x <- x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#View traceplots for multiple parameters
#traceplot(fit.m2, pars = c("Beta0", "sigma_p"))

#Histogram of parameters
posterior <- as.array(fit.m4)
dimnames(posterior) #To see variable names
color_scheme_set("green")
mcmc_hist(posterior, pars = c("Beta0[1]", "Beta0[2]", "Beta0[3]", "Beta0[4]"))

#OBSERVED DATA
obs_data <- coverpercent %>% 
  #filter(Species %in% c("green_algae", "microcoleus", "anabaena_cylindrospermum")) %>% #filter for only 2 sp for now
  filter(year %in% "2022") %>% 
  group_by(field_date, Species) %>% 
  dplyr::summarise(obs_mean = mean(Abundance)) %>% 
  ungroup() %>% 
  mutate(week = rep(seq(1, 13, 2), each = length(unique(Species))))

#Manually calculate mean posteriors for species $ cover, as well as confidence interval
params1 <- as.data.frame(rstan::extract(fit.m2, permuted=FALSE)) %>% 
  select(-c(1:`chain:3.Beta1[4]`)) %>% 
  select(-c(`chain:1.lp__`:`chain:3.lp__`)) %>% 
  mutate(across(1:`chain:3.n[13,4]`, exp)) %>% 
  t 
params2 <- as.data.frame(params1) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(mean = mean(c_across(starts_with("V")), na.rm = TRUE),
                   se_mean = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  mutate(Species = case_when(grepl(",1]", group) ~ 'green_algae',
                             grepl(",2]", group) ~ 'microcoleus',
                             grepl(",3]", group) ~ 'anabaena_cylindrospermum',
                             grepl(",4]", group) ~ 'other_nfixers',
                             grepl("b", group) ~ 'bare_biofilm')) %>% 
  mutate(week = as.numeric(str_extract(group, "[0-9]+")))

#MODEL INCLUDING ALL YEARS--------------------------------------------------------------
#OBSERVED DATA
obs_data_all <- coverpercent %>% 
  group_by(field_date, year, Species) %>% 
  dplyr::summarise(obs_mean = mean(Abundance), obs_SE = calcSE(Abundance)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  mutate(real_week = week(field_date), week = real_week - first(real_week) + 1,
         model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7)) %>% 
  dplyr::filter(Species != "bare_biofilm")

#Manually calculate mean posteriors for species $ cover, as well as confidence interval
params1_all <- as.data.frame(rstan::extract(fit.m4, permuted=FALSE)) %>% 
  dplyr::select(-c(1:`chain:3.Beta_off[4,4]`)) %>% 
  dplyr::select(-c(`chain:1.lp__`:`chain:3.lp__`)) %>% 
  dplyr::select(-c(`chain:1.Ntheta[1]`:`chain:3.Beta[4,4]`)) %>% 
  mutate(across(1:`chain:3.n[4,45]`, exp)) %>% 
  t 

#Set up dataframe to extract week/year info from
yearweek <- alltaxatime %>% 
   pivot_longer(cols = c(green_algae, microcoleus, anabaena_cylindrospermum,
                         bare_biofilm, other_nfixers),
                names_to = "Species", values_to = "mean") %>% 
   mutate(time = rep(seq(45), each = length(unique(Species))))
 
  
params2_all <- as.data.frame(params1_all) %>% 
  rownames_to_column(var="ID") %>% 
  tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
  dplyr::select(-chain) %>% 
  group_by(group) %>% 
  dplyr::summarise(mean = mean(c_across(starts_with("V")), na.rm = TRUE),
                   se_mean = calcSE(c_across(starts_with("V"))),
                   CIlower = quantile(c_across(starts_with("V")), probs = 0.025),
                   CIupper = quantile(c_across(starts_with("V")), probs = 0.975)) %>% 
  mutate(Species = case_when(grepl("1,", group) ~ 'green_algae',
                             grepl("2,", group) ~ 'microcoleus',
                             grepl("3,", group) ~ 'anabaena_cylindrospermum',
                             grepl("4,", group) ~ 'other_nfixers',
                             grepl("b", group) ~ 'bare_biofilm')) %>% 
  mutate(time = as.numeric(ifelse(grepl("b", group), str_extract(group, "[0-9]+"),
         str_extract_all(group, "[0-9]+", simplify = T)[,2]))) %>% 
  left_join(yearweek[,c("uniqueID", "Species", "time")], by = c("Species", "time")) %>% 
  relocate(uniqueID) %>% 
  separate(uniqueID, into = c("year", "week"), sep = "_") %>% 
  mutate(week = as.numeric(week), year = as.numeric(year)) %>% 
  ungroup() %>% 
  left_join(obs_data_all[,c("year", "week", "Species", "real_week")], by = c("year", "week", "Species")) %>% 
  arrange(time) %>% 
  mutate(real_week = ifelse(is.na(real_week), zoo::na.locf(real_week)+1, real_week)) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                              (real_week - 1) * 7 - 1, "week", week_start = 7)) %>% 
  dplyr::filter(Species != "bare_biofilm")




 #FIGURES--------------------------------------------------------------------------------

# #Single year
# ggplot(params2, aes(x = week, y = mean)) + 
#   geom_point(aes(colour = Species), size = 4) + 
#   geom_line(aes(colour = Species), size = 2, alpha = .7) +
#   geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, 
#                   fill = Species), alpha = 0.3) +
#   #geom_errorbar(aes(ymin=mean-se_mean, ymax=mean+se_mean), width=.1) + 
#   geom_point(data = obs_data, aes(x = week, y = obs_mean, shape = Species), size = 2.5) +
#   geom_line(data = obs_data, aes(x = week, y = obs_mean, group = Species),
#             size = .5) +
#   scale_x_continuous(breaks=c(seq(1,13,2))) +
#   labs(x = "Week", y = "Percent Cover (%)", title = "Modeled vs. Observed Abundances, 2022") +
#   labs(color = "Modeled Species", fill = "Modeled Species", shape = "Observed Species")
  
#All years
ggplot(params2_all, aes(x = model_date, y = mean)) + 
  facet_wrap(~year, scales = "free") + 
  geom_point(aes(colour = Species), size = 3) + 
  geom_line(aes(colour = Species), size = 2, alpha = .7) +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, 
                  fill = Species), alpha = 0.3) +
  #geom_errorbar(aes(ymin=mean-se_mean, ymax=mean+se_mean), width=.1) + 
  geom_point(data = obs_data_all, aes(x = model_date, y = obs_mean, shape = Species), size = 2.5) +
  geom_line(data = obs_data_all, aes(x = model_date, y = obs_mean, group = Species),
            size = .5) +
  #scale_x_continuous(breaks=c(seq(1,17,2))) +
  scale_y_continuous(breaks=c(seq(0,100,10))) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Modeled vs. Observed Abundances") +
  labs(color = "Modeled Species", fill = "Modeled Species", shape = "Observed Species")
  
#Pulling out basic numbers
aggregate(mean ~ Species + year, data = params2_all, max)
aggregate(obs_mean ~ Species + year, data = obs_data_all, max)

#Columns plot - need to add barebio back into params2_all for this
ggplot(obs_data_all, aes(x = model_date, y = obs_mean, fill = Species)) +
  facet_wrap(~year, scales = "free") +
  geom_col(position = "fill", width = 5) #+
  #scale_x_continuous(breaks=c(seq(1,17,2))) This was when x = week

#Look at a single year, compare with predictions for that year
subsetallyears <- subset(params2_all, year == 2024)

ggplot(subsetallyears, aes(x = model_date, y = mean)) + 
  geom_point(aes(colour = Species), size = 3) + 
  geom_line(aes(colour = Species), size = 2, alpha = .7) +
  geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`, 
                  fill = Species), alpha = 0.3) +
  #geom_errorbar(aes(ymin=mean-se_mean, ymax=mean+se_mean), width=.1) + 
  geom_point(data = subset(obs_data_all, year == 2024), aes(x = model_date, y = obs_mean, shape = Species), size = 2.5) +
  geom_line(data = subset(obs_data_all, year == 2024), aes(x = model_date, y = obs_mean, group = Species),
            size = .5) +
  scale_y_continuous(breaks=c(seq(0,100,5))) +
  labs(x = "Date", y = "Percent Cover (%)", title = "Modeled vs. Observed Abundances") +
  coord_cartesian(ylim = c(0,35)) +
  labs(color = "Modeled Species", fill = "Modeled Species", shape = "Observed Species")

