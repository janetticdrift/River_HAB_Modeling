###########################
#Data Cleaning: Synthesizing Datasets for Analysis
###########################
#Read in datasets on benthic algal percent cover, microscopy community assembly,
  #nutrient analyses (nitrate, phosphate, and ammonia), temperature,
  #conductivity, discharge, shortwave radiation, and anatoxins. Clean, organize, and 
  #synthesize datasets together into a final analytical framework to be used for 
  #further modeling and analyses.

#Nutrient analyses performed by UNR, temperature and conductivity collected by UNR,
  #discharge rates acquired through USGS National Water Information System (NWIS), 
  #shortwave radiation acquired through NASA North American Land Data Assimilation
  #System (NLDAS), and toxin concentrations measured by SUNY College of Environmental
  #Science and Forestry.

#Packages for cleaning and visualizing data
library(here)
library(tidyverse)
library(dplyr)
library(plyr)
library(ggpubr)
library(devtools)
library(lubridate)
library(patchwork)
library(gridExtra)
library(MASS)
library(slider)

#Package for retrieving USGS discharge data
library(dataRetrieval)

#Packages for retrieving radiation data
library("StreamLightUtils")
library("StreamLight")
library(zoo)
library(ncdf4)
library(CFtime)
library(lattice)
library(httr)

#Read in functions
source(here::here("data_cleaning/Functions.R"))

#Read in river-wide raw data for percent cover by reach and year
percover <- read.csv(here::here("data/percover_byreach.csv")) #2022 and 2023 data
newpercover <- read.csv(here::here("data/SFE ATX % Cover.csv")) #2024 data

#Read in mat community raw data for mat proportion by reach and year
microdata <- read.csv(here::here("data/Target Microscopy.csv"))

#Read in anatoxin data for 2024
atx2223 <- read.csv(here::here("data/cyano_atx.csv"))
atx2024 <- read.csv(here::here("data/cyano_atx_24.csv"))

#Clean new percent cover data to match previous years' formatting
percoverMiranda <- percover %>% 
  dplyr::filter(site == "SFE-M") %>% 
  dplyr::select(!(10:14)) %>%
  dplyr::mutate(field_date = as.Date(field_date)) %>% 
  dplyr::mutate(year = year(field_date)) 

percoverpivot <- percoverMiranda %>% 
  pivot_longer(green_algae:other_nfixers, names_to = "Species", values_to = "Abundance")

newpercoverMiranda <- newpercover %>% 
  dplyr::mutate(site_reach = case_when(Site == "Eel-4S" ~ "SFE-M-1S",
                                Site == "Eel-BUG" ~ "SFE-M-2",
                                Site == "Eel-3UP" ~ "SFE-M-3",
                                Site == "Eel-2UP" ~ "SFE-M-4")) %>% 
  dplyr::rename(green_algae = "GA", microcoleus = "M", anabaena_cylindrospermum = "A",
                other_nfixers = "O", bare_biofilm = "B") %>% 
  dplyr::mutate(field_date = as.Date(format(mdy(Date), '%Y-%m-%d'))) %>% 
  group_by(field_date, site_reach) %>% 
  dplyr::summarise(green_algae = mean(green_algae), microcoleus = mean(microcoleus),
                   anabaena_cylindrospermum = mean(anabaena_cylindrospermum), 
                   other_nfixers = mean(other_nfixers), bare_biofilm = mean(bare_biofilm)) %>% 
  tidyr::separate(site_reach, into = c("site", "M", "reach"), sep="-", remove = FALSE) %>% 
  tidyr::unite("site", c("site", "M"), sep = "-") %>% 
  dplyr::mutate(year = year(field_date)) 

newpercoverpivot <- newpercoverMiranda %>% 
  pivot_longer(green_algae:bare_biofilm, names_to = "Species", values_to = "Abundance")


#Create dataframe that will have Time Step Number and Week Number columns added next 
  #within this file
percoverdata <- rbind(percoverMiranda, newpercoverMiranda)

#Save dataframe of observed percent cover data for plotting modeled vs observed value plots
obs_river_wide_viz <- rbind(percoverpivot, newpercoverpivot)


#Clean microscopy data to tidy slide replicates and consolidate rare species
microscopy_non_avg <- microdata %>% 
  dplyr::filter(grepl("SFE", site_reach)) %>%  #Keep sites that include string "SFE"
  #Updated ID changes: Merge leptolyngbya and geitlerinema
  dplyr::mutate(leptolyngbya_and_geitlerinema = if_else(is.na(leptolyngbya),
                                                        geitlerinema,
                                                        leptolyngbya + geitlerinema)) %>% 
  dplyr::select(!c("leptolyngbya", "geitlerinema")) %>% 
  #Updated ID changes: Move Homoeothrix IDs to Calothrix
  dplyr::mutate(calothrix = if_else(is.na(homoeothrix),
                                    calothrix,
                                    calothrix + homoeothrix)) %>% 
  dplyr::select(!"homoeothrix") %>% 
  #Updated ID changes: Rivularia may be Gloeotrichia pre akinete formation
  dplyr::rename(rivularia_or_early_stage_gloeotrichia = rivularia) %>% 
  #Updated ID changes: Change Lyngbya to Miscellaneous Oscillatoriales
  dplyr::rename(miscellaneous_oscillatoriales = lyngbya) %>% 
  separate(site_reach, into=c("site", "location", "reach"), 
           sep="-") %>% #Split location columns into separate categories
  dplyr::filter(grepl("M", location)) %>% #Keep Miranda sites, remove Standish-Hicky (SH) sites
  dplyr::mutate(field_date = as.Date(field_date, format="%m/%d/%y")) %>%  #Fix date data type
  dplyr::mutate(year = year(field_date)) %>%  #create year column
  relocate(year, .after = field_date) %>%   #reorganize column order
  dplyr::select(!non_algal) %>% #Remove column measuring sediment amount
  dplyr::filter(!reach == 2) %>% #Remove the added reach
  pivot_longer(microcoleus:leptolyngbya_and_geitlerinema, names_to = "Species", values_to = "Abundance")

#Average together slide replicates  
averaged_slides <- microscopy_non_avg %>% 
  dplyr::filter(!slide_rep == "Final") %>% #Do not re-average finalized slide entries
  group_by(field_date, year, site, location, reach, sample_type, date_analyzed, 
           method, Species) %>% 
  dplyr::summarise(Abundance = mean(Abundance)) %>% 
  dplyr::mutate(slide_rep = "Final") %>% 
  relocate(slide_rep, .after = sample_type)

#Isolate already processed (averaged) slides
processed_slides <- microscopy_non_avg %>% 
  dplyr::filter(slide_rep == "Final")

#Bind together averaged replicate slides dataframes into final set
microscopy1 <- rbind(averaged_slides, processed_slides)

#Collapse rare species into one column
#Rare is categorized as any species that has zero occurences in a sample year, for either TM or TAC
non_occurences <- microscopy1 %>% 
  group_by(year, sample_type, Species) %>% 
  dplyr::summarise(Sum = sum(Abundance))

rownum <- which(non_occurences$Sum <= 0) 

rare_species <- non_occurences %>% 
  ungroup() %>% 
  dplyr::slice(rownum)

rare_names <- unique(rare_species$Species)  #Names of the taxa that are rare
rare_names <- c(rare_names, "oscillatoria") #Add taxa that did not appear in 2024

microscopy <- microscopy1 %>% 
  ungroup() %>% 
  pivot_wider(names_from = Species, values_from = Abundance) %>% 
  dplyr::mutate(rare = rowSums(dplyr::select(., rare_names))) %>% 
  dplyr::select(!all_of(rare_names)) %>% 
  dplyr::rename(Anabaena = anabaena_and_cylindrospermum, 
                'Epithemia Diatoms' = e_diatoms,
                Geitlerinema = leptolyngbya_and_geitlerinema,
                'Green Algae' = green_algae, 
                Microcoleus = microcoleus,
                'Non-Epithemia Diatoms' = non_e_diatoms,
                Nostoc = nostoc,
                'Other Coccoids' = other_coccoids,
                Rare = rare) %>% 
  pivot_longer(cols = c("Anabaena":"Rare"), names_to = "Species", values_to = "Abundance")


#Index the dates by timesteps (which timestep in the series) and week numbers 
  #(what number week of the year it is)

#split dataset up into each year
cover_indexdate <- percoverdata %>% 
  group_split(year)

micro_indexdate <- microscopy%>% 
  group_split(year)

year1cover <- cover_indexdate[[1]]
year2cover <- cover_indexdate[[2]]
year3cover <- cover_indexdate[[3]]

year1micro <- micro_indexdate[[1]]
year2micro <- micro_indexdate[[2]]
year3micro <- micro_indexdate[[3]]

#Add time step and week columns
#year 2022
year1_indexdate <- year1cover %>% 
  dplyr::mutate(timestep = dense_rank(field_date)) %>% 
  dplyr::mutate(week = rep(seq(1, 13, 2), times = length(unique(reach))))

year1_indexmicro <- year1micro %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - min(real_week) + 1) %>%
  arrange(reach, field_date) %>% 
  relocate(week, .after = field_date) %>% 
  dplyr::select(!real_week)

#year 2023
year2_indexdate <- year2cover %>% 
  dplyr::mutate(timestep = dense_rank(field_date)) %>% 
  dplyr::mutate(week = timestep)

year2_indexmicro <- year2micro %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - min(real_week) + 1) %>%
  arrange(reach, field_date) %>% 
  relocate(week, .after = field_date) %>% 
  dplyr::select(!real_week)

#year 2024
year3_indexdate <- year3cover %>% 
  dplyr::mutate(timestep = dense_rank(field_date)) %>% 
  arrange(reach) %>% 
  dplyr::mutate(week = rep(seq(1, 17, 2), times = length(unique(reach))))

year3_indexmicro <- year3micro %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - min(real_week) + 1) %>%
  arrange(reach, field_date) %>% 
  relocate(week, .after = field_date) %>% 
  dplyr::select(!real_week)

#River-Wide percent cover binding
cover_indexweek <- rbind(year1_indexdate, year2_indexdate, year3_indexdate)

#Within-Mat community composition binding
micro_indexweek <- rbind(year1_indexmicro, year2_indexmicro, year3_indexmicro) %>% 
  unite(site_reach, c(site, reach), sep = "-", remove =F) %>% 
  arrange(field_date)






#############################################################################################
#Tidy water chemistry data
nut_dat <- read.csv(here::here("data/water_chemistry.csv")) #All years included

#If there is a replacement entry for nitrate or ammonia (usually the minimum 
#detection level), use it instead
nut_data <- nut_dat %>% 
  dplyr::mutate(nitrate_mg_N_L = replace_when(nitrate_mg_N_L, !is.na(nitrate_replace) ~ nitrate_replace),
                ammonia = replace_when(ammonia, !is.na(ammonia_replace) ~ ammonia_replace)) 
                #For any instance where the replacement value is not an NA, replace it with replacement value in the real column

#Subset 2024 measurements and calculate NH4 (ammonium) using ammonia, pH, and temperature 
x <- calculate_NH4(nut_data)[138:173,] %>% 
  dplyr::select(!c(pKa, f))
#Fill back 2024 ammonium calculations into the full dataset
nut_data1 <- nut_data %>% 
  dplyr::slice(1:137)
nut_data2 <- rbind(nut_data1, x)

#Take out an outlier (recording mistake)
nut_data2[85, "nitrate_mg_N_L"] <- NA
#Fix glitch reading from sensor with lowest HOBO estimate
nut_data2[145, "cond_uS_cm"] <- 237 
#Add Dissolved Inorganic Nitrates calculations (sum of nitrate + ammonium)
nut_data3 <- nut_data2 %>% 
  dplyr::mutate(DIN = nitrate_mg_N_L + ammonium_mg_N_L) %>% 
  dplyr::relocate(DIN, .after = ammonium_mg_N_L)
  


#Pull out environmental variables of interest. This cleaned dataset is to be used for 
  #visualizing raw collected variables.
nutrients_raw <- nut_data3 %>% 
  dplyr::filter(site == "SFE-M") %>% 
  dplyr::select(!c("time", "assumed_pH":length(unique(nut_data3)))) %>% #Remove variables that we will not be analyzing
  dplyr::mutate(field_date = as.Date(field_date, format = "%m/%d/%y")) %>% #Update data type
  dplyr::mutate(year = year(field_date)) #Add a grouping column for year
#Average all environmental variables by reach
  #This dataframe is used to create a panel in Figure 1
nutrients_raw_clean <- nutrients_raw %>% 
  group_by(field_date, year) %>% 
  dplyr::summarise(oPhos_ug_P_L = mean(oPhos_ug_P_L), nitrate_mg_N_L = mean(nitrate_mg_N_L, na.rm = T),
                   ammonium_mg_N_L = mean(ammonium_mg_N_L), DIN = mean(DIN, na.rm = T), temp_C = mean(temp_C),
                   cond_uS_cm = mean(cond_uS_cm))
  

#Create similar dataset that interpolates environmental variables for the 
  #weeks in 2022 and 2024 where samples were not collected. Samples in these years 
  #were collected on a biweekly schedule, and 2023 samples were collected weekly

nutrients <- nutrients_raw %>% 
  group_by(year) %>% 
  dplyr::mutate(real_week = week(field_date), week = real_week - first(real_week) + 1) %>% 
  group_by(year) %>% 
  complete(nesting(site_reach, site, reach), week = seq(min(week), max(week), 1L)) %>% #per year week
  ungroup() %>% 
  dplyr::mutate(across(c(temp_C, cond_uS_cm, oPhos_ug_P_L, nitrate_mg_N_L, ammonium_mg_N_L, DIN, real_week), 
                       ~ zoo::na.approx(.x, rule = 2))) %>%  #interpolate env values, and fill in real week NAs
  dplyr::mutate(date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                      (real_week - 1) * 7 - 1, "week", week_start = 7))

#Average environmental variables by reach 
nutrients_avg <- nutrients %>% 
  group_by(date, year) %>% 
  dplyr::summarise(oPhos_ug_P_L = mean(oPhos_ug_P_L), nitrate_mg_N_L = mean(nitrate_mg_N_L),
                   ammonium_mg_N_L = mean(ammonium_mg_N_L), DIN = mean(DIN), temp_C = mean(temp_C),
                   cond_uS_cm = mean(cond_uS_cm))

#Calculate the averae water temperature per year
peryeartempavg <- nutrients_avg %>% 
  group_by(year) %>% 
  dplyr::summarise(mean = mean(temp_C), SE = calcSE(temp_C))  

#Test for normality, and transform data as needed. Parameters that tested normal were not transformed
#Nitrate
shapiro.test(nutrients$nitrate_mg_N_L) #Test for Normality
nutrients$nitrate_mg_N_L <- log(nutrients$nitrate_mg_N_L)

#Phosphate
shapiro.test(nutrients$oPhos_ug_P_L) #Test for Normality
nutrients$oPhos_ug_P_L <- log(nutrients$oPhos_ug_P_L)

#Ammonium
shapiro.test(nutrients$ammonium_mg_N_L) #Test for Normality
nutrients$ammonium_mg_N_L <- log(nutrients$ammonium_mg_N_L)

#DIN
shapiro.test(nutrients$DIN) #Test for Normality
nutrients$DIN <- log(nutrients$DIN)

#Standardize all nutrients to appear on the same scale
stand_nut <- nutrients %>% 
  dplyr::mutate(across(c(oPhos_ug_P_L, nitrate_mg_N_L, ammonium_mg_N_L, DIN, temp_C, cond_uS_cm), 
                ~ scale(.x))) %>% 
  group_by(date, year) %>% 
  dplyr::summarise(oPhos_ug_P_L = mean(oPhos_ug_P_L), nitrate_mg_N_L = mean(nitrate_mg_N_L),
                   ammonium_mg_N_L = mean(ammonium_mg_N_L), DIN = mean(DIN), temp_C = mean(temp_C),
                   cond_uS_cm = mean(cond_uS_cm))

#Visualizing Patterns in Raw Env Data:

#Nitrate
  #Grouped by reach
ggplot(nutrients, aes(x = date, y = nitrate_mg_N_L, colour = reach)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_bw()

ggplot(nutrients_raw_clean, aes(x = field_date, y = nitrate_mg_N_L)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  labs(title = "Nitrate") +
  theme_bw()

#Phosphate
ggplot(nutrients, aes(x = date, y = oPhos_ug_P_L, colour = reach)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  theme_bw()

phosplot <- ggplot(nutrients_avg, aes(x = date, y = oPhos_ug_P_L)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  labs(title = "Phosphate") +
  theme_bw()

#Ammonium
ggplot(nutrients, aes(x = date, y = ammonium_mg_N_L, colour = reach)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  theme_bw()

ammplot <- ggplot(nutrients_avg, aes(x = date, y = ammonium_mg_N_L)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  labs(title = "Ammonium") +
  theme_bw()

#Conductivity
ggplot(nutrients, aes(x = date, y = cond_uS_cm, colour = reach)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  theme_bw()

condplot <- ggplot(nutrients_avg, aes(x = date, y = cond_uS_cm)) +
  facet_wrap(~year, scales = "free_x") +
  geom_point() +
  geom_line() +
  viridis::scale_color_viridis(discrete=TRUE, option="viridis") +
  labs(title = "Conductivity") +
  theme_bw()

cor(nutrients_avg$ammonium_mg_N_L[-c(14:15, 29:30)], anatoxin_data$ATX_all_ug_g)

#Discharge data ---------------------------------------------------------------------
#Acquire discharge flow rates from USGS NWIS

#2022 - startDate = "2022-06-26", endDate = "2022-09-18"
miranda2022 <- renameNWISColumns(readNWISuv(
  siteNumbers = "11476500",
  parameterCd = "00060", #discharge code, cubic feet per second
  startDate = "2022-06-20", #"2021-11-07", for visualizing water year discharge patterns
  endDate = "2022-09-18")) %>% 
  dplyr::mutate(date = as.Date(dateTime)) %>% 
  group_by(date) %>% 
  dplyr::summarise(discharge = mean(Flow_Inst)) %>% 
  dplyr::mutate(discharge = slide_dbl(
    discharge, mean,
    .before = 6, # include previous 6 days into mean for 7 week rolling avg
    .complete = TRUE # calculate mean on full 7-day window
  )) %>%
  dplyr::filter(!is.na(discharge)) %>%
  dplyr::filter(row_number() %% 7 == 1)

#2023 - startDate = "2023-06-20", endDate = "2023-09-24"
miranda2023 <- renameNWISColumns(readNWISuv(
  siteNumbers = "11476500",
  parameterCd = "00060", #discharge code
  startDate = "2023-06-14", #"2022-11-07",
  endDate = "2023-09-25")) %>% 
  dplyr::mutate(date = as.Date(dateTime)) %>% 
  group_by(date) %>% 
  dplyr::summarise(discharge = mean(Flow_Inst)) %>% 
  dplyr::mutate(discharge = slide_dbl(
    discharge, mean,
    .before = 6, # include previous 6 days into mean for 7 week rolling avg
    .complete = TRUE # calculate mean on full 7-day window
  )) %>%
  dplyr::filter(!is.na(discharge)) %>%
  dplyr::filter(row_number() %% 7 == 1)

#2024 - startDate = "2024-06-19", endDate = "2024-10-10"
miranda2024 <- renameNWISColumns(readNWISuv(
  siteNumbers = "11476500",
  parameterCd = "00060", #discharge code
  startDate = "2024-06-13", #"2023-11-02",
  endDate = "2024-10-10")) %>% 
  dplyr::mutate(date = as.Date(dateTime)) %>% 
  group_by(date) %>% 
  dplyr::summarise(discharge = mean(Flow_Inst)) %>% 
  dplyr::mutate(discharge = slide_dbl(
      discharge, mean,
      .before = 6, # include previous 6 days into mean for 7 day avg
      .complete = TRUE # calculate mean on full 7-day window
    )) %>%
  dplyr::filter(!is.na(discharge)) %>%
  dplyr::filter(row_number() %% 7 == 1)

discharge <- rbind(miranda2022, miranda2023, miranda2024) %>% 
  dplyr::mutate(year = factor(year(date))) %>% 
  # dplyr::mutate(season_year = if_else(month(date) >= 11, year(date) + 1, year(date)), #Make Nov and Dec of previous year belong to the same season year, if visualizing water year 
  #               season_year = factor(season_year)) %>%
  dplyr::mutate(fake_date = make_date(year = min(year(date)), day = day(date), month = month(date)),
                fake_date = if_else(month(date) >= 11, fake_date - years(1), fake_date)) %>% 
  dplyr::mutate(log_discharge = log(discharge)) %>% #log-transform
  dplyr::mutate(stand_discharge = c(scale(log_discharge))) 
  #dplyr::mutate(year = as.numeric(as.character(year))) #Use when using season_year


#Visualization of raw discharge patterns
ggplot(discharge, aes(x = fake_date, y = discharge, group = year, color = year)) +
  geom_point() +
  geom_line(size = 1) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b")+ #b = month?
  labs(x = "Date")

#Zoom in on the dates past the spring peak
ggplot(discharge, aes(x = fake_date, y = log_discharge, color = year)) +
  geom_point() +
  geom_line() +
  coord_cartesian(xlim = as.Date(c('2022-06-01', '2022-11-01')), ylim = c(0,7))

#############################################################################################
#Import and tidy photosynthetically active radiation (PAR) data

#The commented out code is no longer UTD as of November 17th 2025, as the data source is 
  #housed in a new repository now
#I downloaded the PAR data within relevant dates from last load, and now archived on GitHub instead

# source("/Users/jld/Documents/Github/River_HAB_Modeling/data_cleaning/Functions.R")
# 
# #Process and format NLDAS data for last two months of 2024. SW = shortwaves
# NLDAS_sw <- get_NLDASv20_datarod(
#   start_date = "2022-06-26",
#   end_date = "2024-10-10",
#   lat = 40.198173,
#   lon = -123.775930,
#   var = "SWdown"
# )
# 
# #Separate data out into the dates used per year
# PAR <- NLDAS_sw %>% 
#   dplyr::rename(radiation = value) %>% #Metric is SW_W_m_2
#   separate(datetime, c("date", "time"), sep = " ") %>% 
#   dplyr::group_by(date) %>% 
#   dplyr::summarise(radiation = mean(radiation))

PAR <- read.csv(here::here("data/PAR.csv"))

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
swradiation_raw <- rbind(PAR2022, PAR2023, PAR2024) %>% 
  dplyr::mutate(year = factor(year(date))) %>% 
  dplyr::mutate(fake_date = make_date(year = min(year(date)), day = day(date), month = month(date))) %>% 
  dplyr::mutate(date = as.Date(date),
                year = as.numeric(as.character(year)))

#Scale radiation data
swradiation <- swradiation_raw %>% 
  dplyr::mutate(stand_rad = c(scale(radiation)))

#Visualization of raw radiation data
ggplot(swradiation, aes(x = fake_date, y = radiation, group = year, color = year)) +
  geom_point() +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") #b = month?

#############################################################################################
                        #Plot Nitrate and Discharge Together for Figure 2

#Merge together nitrate and discharge data
nitrate.discharge <- cbind(nutrients_raw_clean[, c("year", "nitrate_mg_N_L")], 
                           discharge[, c("discharge", "fake_date"), drop = FALSE]) %>% 
  dplyr::rename(nitrate = nitrate_mg_N_L) %>% 
  #To make the "fake" dates within the same year match perfectly with field dates, replace the last date 
    #(10-09-2022) with the real field date (10-13-2022)
  dplyr::mutate(fake_date = replace(fake_date, 45, "2022-10-13")) #45 means at index position 45

scale_factor <- 5000 #Rough estimate to scale up nitrate by 

envplot <- ggplot(nitrate.discharge, aes(x = fake_date)) +
  facet_wrap(~year) +
  geom_line(aes(y = discharge, color = "Discharge"), size = 1) +
  geom_point(aes(y = discharge, color = "Discharge")) +
  geom_line(aes(y = nitrate*scale_factor, color = "Nitrate"), size = 1) +
  geom_point(aes(y = nitrate*scale_factor, color = "Nitrate")) +
  scale_y_continuous(
    name = "Discharge (cfs)",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Nitrate (mg N/L)"
    )
  ) +
  scale_color_manual(name = "Env. Variable", 
                     values = c("Discharge" = "#813B9A",
                                "Nitrate" = "#1a7531")) +
  labs(x = "Date") +
  theme_bw()

#############################################################################################
                                #Tidy Anatxoin Concentration data.

#Clean and subset anatoxin data from 2022 and 2023
atx2223clean <- atx2223 %>% 
  dplyr::filter(grepl("SFE", site_reach)) %>%  #Keep sites that include string "SFE" in site col
  dplyr::select(!c(site, site_reach)) %>% 
  dplyr::select(!c(Chla_ug_g:12)) %>%  #Remove toxins that weren't analyzed in 2024
  dplyr::mutate(field_date = as.Date(field_date))

#Clean and subset anatoxin data from 2024
atx24clean <- atx2024 %>% 
  dplyr::filter(grepl("SFE", site)) %>%  #Keep sites that include string "SFE" in site col
  dplyr::mutate(is_dup = grepl("Duplicate", Sample)) %>% #create empty col that stores duplicate info
  dplyr::mutate(across(.cols = where(is.numeric), #only target numeric columns
                .fns = ~ if_else(is_dup, (. + lag(.)) / 2, .))) %>%  #.row + preceding .row / 2. else, keep row same
  dplyr:: mutate(across(where(is.numeric),
                ~ if_else(replace_na(lead(is_dup), FALSE), lead(.), .))) %>% #if next row has is_dup=T, replace current row with next row's values. 
                                                                             #replace_NA says to NOT replace rows with NA, since the last row does not have a next row for lead() to work on it returns NAs
  dplyr::filter(!is_dup) %>% #remove duplicate rows where is_dup = T
  dplyr::filter(!grepl("Var Reps", Sample)) %>% #Remove extra analysis replicates
  dplyr::select(!c(is_dup, Sample, site, ESF)) %>%  #remove duplicate ID col and Sample col, and unnecessary ESF sample ID column
  dplyr::select(!c(Total_ATXs:Det_Limits_MCs, dhHTXa_ug_g)) %>%   #remove toxins that weren't analyzed in '22 or '23
  dplyr::mutate(field_date = as.Date(format(mdy(field_date), '%Y-%m-%d'))) %>% 
  arrange(field_date)
  
#combine 2024 data with 2022,2023

toxindf <- rbind(atx2223clean, atx24clean) %>% 
  pivot_longer(4:7, names_to = "anatoxins", values_to = "concentration") %>% 
  group_by(field_date, reach, sample_type, anatoxins) %>% 
  dplyr::summarise(concentration = mean(concentration)) %>%  #For reaches with multiple samples, take the average
  dplyr::mutate(year = year(field_date)) %>% 
  pivot_wider(names_from = "anatoxins", values_from = "concentration") %>% 
  dplyr::mutate(sample_type = case_when(sample_type=="TM" ~ "Microcoleus",
                                        sample_type=="TAC" ~ "Anabaena"))
