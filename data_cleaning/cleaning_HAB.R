#Data cleaning of SFE algae cover data per year

library(plyr)
library(tidyverse)
library(ggpubr)
library(gjam)
library(devtools)
library(here)
library(lubridate)
library("StreamLightUtils")
library("StreamLight")

percover1 <- read.csv(here::here("data/percover_byreach.csv")) #2022 and 2023 data
newpercover <- read.csv(here::here("data/SFE ATX % Cover.csv")) #2024 data

#Clean new data to match previous year formatting

percover <- percover1 %>% 
  filter(site == "SFE-M") %>% 
  select(!(10:14)) %>%
  mutate(field_date = as.Date(field_date)) %>% 
  mutate(year = year(field_date))

percoverplot <- percover %>% 
  pivot_longer(green_algae:other_nfixers, names_to = "Species", values_to = "Abundance")

cleanpercover <- newpercover %>% 
  mutate(site_reach = case_when(Site == "Eel-4S" ~ "SFE-M-1S",
                                Site == "Eel-BUG" ~ "SFE-M-2",
                                Site == "Eel-3UP" ~ "SFE-M-3",
                                Site == "Eel-2UP" ~ "SFE-M-4")) %>% 
  dplyr::rename(green_algae = "GA", microcoleus = "M", anabaena_cylindrospermum = "A",
                other_nfixers = "O", bare_biofilm = "B") %>% 
  mutate(field_date = as.Date(format(mdy(Date), '%Y-%m-%d'))) %>% 
  group_by(field_date, site_reach) %>% 
  dplyr::summarise(green_algae = mean(green_algae), microcoleus = mean(microcoleus),
                   anabaena_cylindrospermum = mean(anabaena_cylindrospermum), 
                   other_nfixers = mean(other_nfixers), bare_biofilm = mean(bare_biofilm)) %>% 
  tidyr::separate(site_reach, into = c("site", "M", "reach"), sep="-", remove = FALSE) %>% 
  tidyr::unite("site", c("site", "M"), sep = "-") %>% 
  mutate(year = year(field_date))

cleanpercoverplot <- cleanpercover %>% 
  pivot_longer(green_algae:bare_biofilm, names_to = "Species", values_to = "Abundance")

allcoverdataplot <- rbind(percoverplot, cleanpercoverplot)
allcoverdata <- rbind(percover, cleanpercover)

#############################################################################################
#Index date by timesteps: 1, 2, 3... n

#split dataset up into each year
cover_indexdate <- allcoverdata %>% 
  group_split(year)

year1cover <- cover_indexdate[[1]]
year2cover <- cover_indexdate[[2]]
year3cover <- cover_indexdate[[3]]

#year 2022
year1_indexdate <- year1cover %>% 
  mutate(timestep = dense_rank(field_date)) %>% 
  mutate(week = rep(seq(1, 13, 2), times = length(unique(reach))))

#year 2023
year2_indexdate <- year2cover %>% 
  mutate(timestep = dense_rank(field_date)) %>% 
  mutate(week = timestep)

#year 2024
year3_indexdate <- year3cover %>% 
  mutate(timestep = dense_rank(field_date)) %>% 
  arrange(reach) %>% 
  mutate(week = rep(seq(1, 17, 2), times = length(unique(reach))))

cover_indexweek <- rbind(year1_indexdate, year2_indexdate, year3_indexdate)

#############################################################################################
#Tidy water chemistry data
nut_data <- read.csv(here::here("data/water_chemistry.csv")) #All years included

#Function for calculating ammonium from ammonia
calculate_NH4 <- function(df) {
  df <- df %>%
    dplyr::mutate(
      pKa = 0.09018 + 2727.92 / (temp_C + 273.15),
      f = 1 / (10^(pKa - pH) + 1),
      ammonium_mg_N_L = round((1 - f) * ammonia, 5)
    )
}
#Subset out the last year
x <- calculate_NH4(nut_data)[138:173,] %>% 
  select(!c(pKa, f))
#Fill back in the full dataset
nut_data <- nut_data %>% 
  slice(1:137)
nut_data <- rbind(nut_data, x)
nut_data[85, "nitrate_mg_N_L"] <- NA #Take out an outlier



#Pull out variables of interest
nutrients <- nut_data %>% 
  filter(site == "SFE-M") %>% 
  select(!c("time", 15:length(unique(nut_data)))) %>% #remove uninteresting nutrients
  mutate(field_date = as.Date(field_date, format = "%m/%d/%y")) %>% 
  mutate(year = year(field_date)) %>% 
  group_by(year) %>% 
  mutate(real_week = week(field_date), week = real_week - first(real_week) + 1) %>% 
  group_by(year) %>% 
  complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% #per year week
  ungroup() %>% 
  mutate(across(c(temp_C, cond_uS_cm, oPhos_ug_P_L, nitrate_mg_N_L, ammonium_mg_N_L, real_week), 
                ~ zoo::na.approx(.x, rule = 2))) %>%  #interpolate env values, and fill in real week NAs
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))

stand_nut <- nutrients %>% 
  mutate(across(c(oPhos_ug_P_L, nitrate_mg_N_L, ammonium_mg_N_L, temp_C, cond_uS_cm), 
                ~ scale(.x))) %>% 
  group_by(model_date, year) %>% 
  dplyr::summarise(oPhos_ug_P_L = mean(oPhos_ug_P_L), nitrate_mg_N_L = mean(nitrate_mg_N_L),
                   ammonium_mg_N_L = mean(ammonium_mg_N_L),temp_C = mean(temp_C),
                   cond_uS_cm = mean(cond_uS_cm))


#Averaged by reach 
nutrients_avg <- nutrients %>% 
  group_by(model_date, year) %>% 
  dplyr::summarise(oPhos_ug_P_L = mean(oPhos_ug_P_L), nitrate_mg_N_L = mean(nitrate_mg_N_L),
                   ammonium_mg_N_L = mean(ammonium_mg_N_L), temp_C = mean(temp_C),
                   cond_uS_cm = mean(cond_uS_cm))

#Prelim plot: Nitrate
ggplot(nutrients, aes(x = model_date, y = nitrate_mg_N_L, colour = reach)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_bw()

ggplot(nutrients_avg, aes(x = model_date, y = nitrate_mg_N_L)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_bw()

#Prelim plot: Phosphate
ggplot(nutrients, aes(x = model_date, y = oPhos_ug_P_L, colour = reach)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  theme_bw()

ggplot(nutrients_avg, aes(x = model_date, y = oPhos_ug_P_L)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  theme_bw()

#Prelim plot: Ammonium
ggplot(nutrients, aes(x = model_date, y = ammonium_mg_N_L, colour = reach)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  theme_bw()

ggplot(nutrients_avg, aes(x = model_date, y = ammonium_mg_N_L)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  theme_bw()

#Prelim plot: Conductivity
ggplot(nutrients, aes(x = model_date, y = cond_uS_cm, colour = reach)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  theme_bw()

ggplot(nutrients_avg, aes(x = model_date, y = cond_uS_cm)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  theme_bw()


#Discharge data ---------------------------------------------------------------------

#2022 - startDate = "2022-06-26", endDate = "2022-09-18"
miranda2022 <- renameNWISColumns(readNWISuv(
  siteNumbers = "11476500",
  parameterCd = "00060", #discharge code, cubic feet per second!
  startDate = "2022-06-26",
  endDate = "2022-09-18")) %>% 
  mutate(date = as.Date(dateTime)) %>% 
  group_by(date) %>% 
  dplyr::summarise(discharge = mean(Flow_Inst)) %>% 
  filter(row_number() %% 7 == 1)

#2023 - startDate = "2023-06-20", endDate = "2023-09-24" but end it a couple days later
miranda2023 <- renameNWISColumns(readNWISuv(
  siteNumbers = "11476500",
  parameterCd = "00060", #discharge code
  startDate = "2023-06-20",
  endDate = "2023-09-25")) %>% 
  mutate(date = as.Date(dateTime)) %>% 
  group_by(date) %>% 
  dplyr::summarise(discharge = mean(Flow_Inst)) %>% 
  filter(row_number() %% 7 == 1)

#2024 - startDate = "2024-06-19", endDate = "2024-10-10"
miranda2024 <- renameNWISColumns(readNWISuv(
  siteNumbers = "11476500",
  parameterCd = "00060", #discharge code
  startDate = "2024-06-19",
  endDate = "2024-10-10")) %>% 
  mutate(date = as.Date(dateTime)) %>% 
  group_by(date) %>% 
  dplyr::summarise(discharge = mean(Flow_Inst)) %>% 
  filter(row_number() %% 7 == 1)

discharge <- rbind(miranda2022, miranda2023, miranda2024) %>% 
  mutate(year = factor(year(date))) %>% 
  mutate(fake_date = make_date(year = min(year(date)), day = day(date), month = month(date))) %>% 
  mutate(log_discharge = log(discharge)) %>% 
  mutate(stand_discharge = c(scale(discharge)))


#Quick plot of discharge data
ggplot(discharge, aes(x = fake_date, y = log_discharge, color = year)) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") #b = month?

#Zoom in on the dates past the spring peak
ggplot(discharge, aes(x = fake_date, y = log_discharge, color = year)) +
  geom_point() +
  geom_line() +
  coord_cartesian(xlim = as.Date(c('2022-06-01', '2022-11-01')), ylim = c(0,7))

#############################################################################################
#Import and tidy photosynthetically active radiation (PAR) data

#Only need to run this download code once
#working_dir <- "/Users/jld/Documents/Github/River_HAB_Modeling/data"
# NLDAS_DL(save_dir = working_dir, Site_ID = "sfkeel_mir", Lat = 40.198173, 
#          Lon = -123.775930, startDate = "2022-06-29")

NLDAS_processed <- NLDAS_proc(working_dir, "sfkeel_mir")[["sfkeel_mir"]]
             
             
