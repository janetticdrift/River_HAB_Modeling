###########################
#Microscopy initial visualizations
###########################

library(tidyverse)
library(ggplot2)
library(here)
library(lubridate)

#Read in 2022 and 2023 data
microdata <- read.csv(here::here("data/Target Microscopy.csv"))

#Clean data for analysis
microscopy <- microdata %>% 
  filter(grepl("SFE", site_reach)) %>%  #Keep sites that include string "SFE"
  separate(site_reach, into=c("site", "location", "reach"), 
           sep="-") %>% #Split location columns into separate categories
  filter(grepl("M", location)) %>% #Keep Miranda sites, remove Steelhead (SH) sites
  mutate(field_date = as.Date(field_date, format="%m/%d/%y")) %>%  #Fix date data type
  mutate(year = year(field_date)) %>%  #create year column
  relocate(year, .after = field_date) %>%   #reorganize column order
  pivot_longer(microcoleus:aphanothece, names_to = "Species", values_to = "Abundance") #Change layout of dataframe

#Plot timeseries of within-mat species per year
# ggplot(microscopy, aes(x = field_date, y = Abundance, fill = Species)) +
#   facet_wrap(~year, scales = "free") +
#   geom_col(position = "fill") Can't see the detail well, subset per year

#2022
ggplot(subset(microscopy, year %in% 2022), aes(x = field_date, y = Abundance, fill = Species)) +
  geom_col(position = "fill") +
  labs(title = "2022") +
  theme(legend.position="bottom")

#2023
ggplot(subset(microscopy, year %in% 2023), aes(x = field_date, y = Abundance, fill = Species)) +
  geom_col(position = "fill") +
  labs(title = "2023")
