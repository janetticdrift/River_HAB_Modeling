#Data cleaning of SFE algae cover data per year

#Packages for cleaning data
library(plyr)
library(tidyverse)
library(ggpubr)
library(gjam)
library(devtools)
library(here)
library(lubridate)
library(dataRetrieval)

#Packages for radiation data
library("StreamLightUtils")
library("StreamLight")
library(zoo)
library(ncdf4)
library(CFtime)
library(lattice)
library(httr)

#Read in river-wide raw data for percent cover by reach and year
percover1 <- read.csv(here::here("data/percover_byreach.csv")) #2022 and 2023 data
newpercover <- read.csv(here::here("data/SFE ATX % Cover.csv")) #2024 data

#Read in mat community raw data for mat proportion by reach and year
microdata <- read.csv(here::here("data/Target Microscopy.csv"))

#Read in anatoxin data for 2024
atx2223 <- read.csv(here::here("data/cyano_atx.csv"))
atx2024 <- read.csv(here::here("data/cyano_atx_24.csv"))

#Clean new percent cover data to match previous year formatting

percover <- percover1 %>% 
  dplyr::filter(site == "SFE-M") %>% 
  dplyr::select(!(10:14)) %>%
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


#Clean microscopy data to tidy slide replicates and consolidate rare species
microscopy_non_avg <- microdata %>% 
  filter(grepl("SFE", site_reach)) %>%  #Keep sites that include string "SFE"
  separate(site_reach, into=c("site", "location", "reach"), 
           sep="-") %>% #Split location columns into separate categories
  filter(grepl("M", location)) %>% #Keep Miranda sites, remove Standish-Hicky (SH) sites
  mutate(field_date = as.Date(field_date, format="%m/%d/%y")) %>%  #Fix date data type
  mutate(year = year(field_date)) %>%  #create year column
  relocate(year, .after = field_date) %>%   #reorganize column order
  select(!non_algal) %>% #Remove column measuring sediment amount
  filter(!reach == 2) %>% #Remove the added reach
  pivot_longer(microcoleus:aphanothece, names_to = "Species", values_to = "Abundance")

#Average together slide replicates  
averaged_slides <- microscopy_non_avg %>% 
  filter(!slide_rep == "Final") %>%
  group_by(field_date, year, site, location, reach, sample_type, date_analyzed, 
           method, Species) %>% 
  dplyr::summarise(Abundance = mean(Abundance)) %>% 
  mutate(slide_rep = "Final") %>% 
  relocate(slide_rep, .after = sample_type)

#Pull out already processed slides
processed_slides <- microscopy_non_avg %>% 
  filter(slide_rep == "Final")

#Bind together dataframes into final set
microscopy1 <- rbind(averaged_slides, processed_slides)

#Collapse rare species into one column
non_occurences <- microscopy1 %>% 
  group_by(year, sample_type, Species) %>% 
  dplyr::summarise(Sum = sum(Abundance))

rownum <- which(non_occurences$Sum <= 0) #Which samples have no occurences per each target and each year?

rare_species <- non_occurences %>% 
  #dplyr::filter(year == "2022") %>% 
  ungroup() %>% 
  dplyr::slice(rownum)

rare_names <- unique(rare_species$Species) #Names of the species that are rare

microscopy <- microscopy1 %>% 
  ungroup() %>% 
  pivot_wider(names_from = Species, values_from = Abundance) %>% 
  dplyr::mutate(rare = rowSums(select(., rare_names))) %>% 
  dplyr::select(!rare_names)

#############################################################################################
#Index date by timesteps and week numbers: 1, 2, 3... n

#split dataset up into each year
cover_indexdate <- allcoverdata %>% 
  group_split(year)

micro_indexdate <- microscopy%>% 
  group_split(year)

year1cover <- cover_indexdate[[1]]
year2cover <- cover_indexdate[[2]]
year3cover <- cover_indexdate[[3]]

year1micro <- micro_indexdate[[1]]
year2micro <- micro_indexdate[[2]]

#Create function for assigning week numbers
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
year1_indexdate <- year1cover %>% 
  mutate(timestep = dense_rank(field_date)) %>% 
  mutate(week = rep(seq(1, 13, 2), times = length(unique(reach))))

year1_indexmicro <- year1micro %>% 
  mutate(timestep = dense_rank(field_date)) %>% 
  arrange(reach, field_date) %>% 
  mutate(week = case_when(
    timestep == 1 ~ 1,
    timestep == 2 ~ 3,
    timestep == 3 ~ 5,
    timestep == 4 ~ 7,
    timestep == 5 ~ 9,
    timestep == 6 ~ 11,
    timestep == 7 ~ 13,
    timestep == 8 ~ 15))

#year 2023
year2_indexdate <- year2cover %>% 
  mutate(timestep = dense_rank(field_date)) %>% 
  mutate(week = timestep)

year2_indexmicro <- year2micro %>% 
  mutate(timestep = dense_rank(field_date)) %>% 
  arrange(reach, field_date) %>% 
  mutate(week = timestep)

#year 2024
year3_indexdate <- year3cover %>% 
  mutate(timestep = dense_rank(field_date)) %>% 
  arrange(reach) %>% 
  mutate(week = rep(seq(1, 17, 2), times = length(unique(reach))))

#River-wide percent cover binding
cover_indexweek <- rbind(year1_indexdate, year2_indexdate, year3_indexdate)

#Mat community proportion binding
micro_indexweek <- rbind(year1_indexmicro, year2_indexmicro) %>% 
  unite(site_reach, c(site, reach), sep = "-", remove =F)






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
  dplyr::select(!c(pKa, f))
#Fill back in the full dataset
nut_data <- nut_data %>% 
  dplyr::slice(1:137)
nut_data <- rbind(nut_data, x)
nut_data[85, "nitrate_mg_N_L"] <- NA #Take out an outlier
nut_data[145, "cond_uS_cm"] <- 237 #Fix glitch reading from sensor with lowest HOBO estimate



#Pull out variables of interest
nutrients <- nut_data %>% 
  dplyr::filter(site == "SFE-M") %>% 
  dplyr::select(!c("time", 15:length(unique(nut_data)))) %>% #remove uninteresting nutrients
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
  dplyr::filter(row_number() %% 7 == 1)

#2023 - startDate = "2023-06-20", endDate = "2023-09-24" but end it a couple days later
miranda2023 <- renameNWISColumns(readNWISuv(
  siteNumbers = "11476500",
  parameterCd = "00060", #discharge code
  startDate = "2023-06-20",
  endDate = "2023-09-25")) %>% 
  mutate(date = as.Date(dateTime)) %>% 
  group_by(date) %>% 
  dplyr::summarise(discharge = mean(Flow_Inst)) %>% 
  dplyr::filter(row_number() %% 7 == 1)

#2024 - startDate = "2024-06-19", endDate = "2024-10-10"
miranda2024 <- renameNWISColumns(readNWISuv(
  siteNumbers = "11476500",
  parameterCd = "00060", #discharge code
  startDate = "2024-06-19",
  endDate = "2024-10-10")) %>% 
  mutate(date = as.Date(dateTime)) %>% 
  group_by(date) %>% 
  dplyr::summarise(discharge = mean(Flow_Inst)) %>% 
  dplyr::filter(row_number() %% 7 == 1)

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

source("/Users/jld/Documents/Github/River_HAB_Modeling/data_cleaning/R Functions/Hydrology Data Rods.R")

#Process and format NLDAS data for last two months of 2024. SW = shortwaves
NLDAS_sw <- get_NLDASv20_datarod(
  start_date = "2022-06-26",
  end_date = "2024-10-10",
  lat = 40.198173,
  lon = -123.775930,
  var = "SWdown"
)

#Separate data out into the dates used per year
PAR <- NLDAS_sw %>% 
  dplyr::rename(radiation = value) %>% #Metric is SW_W_m_2
  separate(datetime, c("date", "time"), sep = " ") %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarise(radiation = mean(radiation))

#Remove non-field season dates for each year
PAR2022 <- PAR %>% 
  dplyr::filter(between(date, "2022-06-26", "2022-09-18")) %>% 
  dplyr::filter(row_number() %% 7 == 1)
PAR2023 <- PAR %>% 
  dplyr::filter(between(date, "2023-06-20", "2023-09-26")) %>% 
  dplyr::filter(row_number() %% 7 == 1)
PAR2024 <- PAR %>% 
  dplyr::filter(between(date, "2024-06-19", "2024-10-10")) %>% 
  dplyr::filter(row_number() %% 7 == 1)

#Bind together yearly PAR data
swradiation <- rbind(PAR2022, PAR2023, PAR2024) %>% 
  mutate(year = factor(year(date))) %>% 
  mutate(fake_date = make_date(year = min(year(date)), day = day(date), month = month(date))) %>% 
  mutate(stand_rad = c(scale(radiation)))

#Quick plot of radiation data
ggplot(swradiation, aes(x = fake_date, y = radiation, color = year)) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") #b = month?


#############################################################################################
#Tidy ATX data.

atx2223clean <- atx2223 %>% 
  filter(grepl("SFE", site_reach)) %>%  #Keep sites that include string "SFE" in site col
  select(!site) %>% 
  select(!c(Chla_ug_g:13)) #Remove toxins that weren't analyzed in 2024

#WHAT ARE THE EXT SAMPLES--I ignored them for now!
atx24clean <- atx2024 %>% 
  filter(grepl("SFE", site_reach)) %>%  #Keep sites that include string "SFE" in site col
  mutate(is_dup = grepl("Duplicate", Sample)) %>% #create col that stores duplicate info
  mutate(across(.cols = where(is.numeric), #only target numeric columns
                .fns = ~ if_else(is_dup, (. + lag(.)) / 2, .))) %>%  #.row + preceding .row / 2. else, keep row same
  mutate(across(where(is.numeric),
                ~ if_else(replace_na(lead(is_dup), FALSE), lead(.), .))) %>% #if next row has is_dup=T, replace current row with next row's values. 
                                                                             #replace_NA says to NOT replace rows with NA, since the last row does not have a next row for lead() to work on it returns NAs
  filter(!is_dup) %>% #remove old duplicate rows
  select(!c(is_dup, Type, Sample)) %>%  #remove duplicate ID col
  dplyr::rename(field_date = Date) %>% 
  mutate(sample_type = "TM") %>% 
  select(!c(Total_ATXs:Det_Limits_MCs, dhHTXa_ug_g)) #remove toxins that weren't analyzed in '22,'23

#combine 2024 data with 2022,2023
atx <- rbind(atx2223clean, atx24clean) %>% 
  pivot_longer(5:8, names_to = "anatoxins", values_to = "concentration") %>% 
  group_by(field_date, reach, sample_type, anatoxins) %>% 
  dplyr::summarise(concentration = mean(concentration)) %>%  #For reaches with multiple samples, average
  mutate(field_date = as.Date(field_date)) %>% 
  mutate(year = year(field_date))

#Plot data
ggplot(subset(atx, sample_type %in% "TM"), aes(x = field_date, y = concentration, color = anatoxins)) +
  facet_grid(reach~year, scales = "free") + #facet_grid for multiple variables
  geom_point() +
  geom_line()
#Row 73 in atx is the 102 concentration of total atx in 2022 reach 4



#----------------------------------------
#Calculation of average temperature
calcSE<-function(x){
  x <- x[is.na(x)==F] 
  sd(x)/sqrt(length(x))
}
calcSE(nutrients_avg$temp_C)
