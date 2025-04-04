###############
#Look at raw data
###############

#Library
library(plyr)
library(tidyverse)
library(ggpubr)
library(gjam)
library(devtools)
library(here)
library(lubridate)

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

#Save data
saveRDS(
  allcoverdataplot, 
  file = here::here("data/allcoverdataplot.rds")
) 


#Plot timeseries of 2024 data
ggplot(cleanpercoverplot, aes(x = field_date, y = Abundance, color = Species)) +
  facet_wrap(vars(reach)) +
  geom_point() +
  geom_line()

#Plot timeseries of all data
ggplot(allcoverdataplot, aes(x = field_date, y = Abundance, color = Species)) +
  facet_wrap(vars(reach)) +
  geom_point() +
  geom_line(aes(group = interaction(Species, year(field_date)))) +
  scale_x_date(guide = guide_axis(n.dodge = 2))

#Zoom into any reach you want
ggplot(subset(allcoverdataplot, reach %in% "1S"), 
       aes(x = field_date, y = Abundance, color = Species)) +
  facet_wrap(vars(reach)) +
  geom_point() +
  geom_line(aes(group = interaction(Species, year(field_date))))

#Proportion table
ggplot(allcoverdataplot, aes(x = field_date, y = Abundance, fill = Species)) +
  facet_wrap(~year, scales = "free") +
  geom_col(position = "fill")



#############################################################################################
#Index date by 1, 2, 3.....

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

#Average across reaches by date
avg_cover <- cover_indexweek %>% 
  group_by(field_date, timestep, week) %>% 
  dplyr::summarise(Avg_greenalgae = mean(green_algae), Avg_microcoleus = mean(microcoleus),
                   Avg_anabaena = mean(anabaena_cylindrospermum), Avg_bare = mean(bare_biofilm),
                   Avg_nfixers = mean(other_nfixers))
  
  
  
  