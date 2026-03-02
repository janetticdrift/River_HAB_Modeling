###########################
#Microscopy initial visualizations
###########################

library(tidyverse)
library(ggplot2)
library(here)
library(lubridate)
library(forcats)
library(RColorBrewer)

#Read in 2022 and 2023 data
microdata <- read.csv(here::here("data/Target Microscopy.csv"))

#Clean data for analysis
microscopy_non_avg <- microdata %>% 
  dplyr::filter(grepl("SFE", site_reach)) %>%  #Keep sites that include string "SFE"
  separate(site_reach, into=c("site", "location", "reach"), 
           sep="-") %>% #Split location columns into separate categories
  dplyr::filter(grepl("M", location)) %>% #Keep Miranda sites, remove Standish-Hicky (SH) sites
  mutate(field_date = as.Date(field_date, format="%m/%d/%y")) %>%  #Fix date data type
  mutate(year = year(field_date)) %>%  #create year column
  relocate(year, .after = field_date) %>%   #reorganize column order
  dplyr::select(!non_algal) %>% #Remove column measuring sediment amount
  pivot_longer(microcoleus:aphanothece, names_to = "Species", values_to = "Abundance")    #Change layout of dataframe
  
  #Average together slide replicates  
  averaged_slides <- microscopy_non_avg %>% 
    dplyr::filter(!slide_rep == "Final") %>%
    group_by(field_date, year, site, location, reach, sample_type, date_analyzed, 
             method, Species) %>% 
    dplyr::summarise(Abundance = mean(Abundance)) %>% 
    mutate(slide_rep = "Final") %>% 
    relocate(slide_rep, .after = sample_type)
  #Pull out already processed slides
  processed_slides <- microscopy_non_avg %>% 
    filter(slide_rep == "Final")
  #Bind together dataframes
  microscopy1 <- rbind(averaged_slides, processed_slides)
  
  
#Do TAC samples all actually contain Anabaena?
TAC <- microscopy1 %>% 
  filter(sample_type == "TAC" & Species == "anabaena_and_cylindrospermum")

which(TAC$Abundance < 10) #Are there any instances where Abundance is small? FALSE = no
which(TAC$Abundance < 20)
  #Some instances where Anabaena is less than 20% in TAC sample

#Are there any species that never appear?
non_occurences <- microscopy1 %>% 
  group_by(year, sample_type, Species) %>% 
  dplyr::summarise(Sum = sum(Abundance))

rownum <- which(non_occurences$Sum <= 0) #Which samples have no occurences per each target and year?

rare_species <- non_occurences %>% 
  ungroup() %>% 
  dplyr::slice(rownum) #Leave in row numbers identified above as having no occurence

rare_names <- unique(rare_species$Species) #Names of the species that are rare
#rare_names <- rare_names[rare_names != "gloeotrichia"]

#Remove non-occuring species, and combine those that are super rare
microscopy <- microscopy1 %>% 
  ungroup() %>% 
  pivot_wider(names_from = Species, values_from = Abundance) %>% 
  dplyr::mutate(rare = rowSums(dplyr::select(., rare_names))) %>% 
  dplyr::select(!rare_names) %>%  #Remove rare columns now that they have been consolidated
  pivot_longer(anabaena_and_cylindrospermum:rare, names_to = "Species", values_to = "Abundance") %>% 
  dplyr::mutate(Species = case_match(Species,
                                     "anabaena_and_cylindrospermum" ~ "Anabaena", 
                                     "e_diatoms" ~ "Epithemia Diatoms", 
                                     "geitlerinema" ~ "Geitlerinema", 
                                     "green_algae" ~ "Green Algae", 
                                     "leptolyngbya" ~ "Leptolyngbya", 
                                     "microcoleus" ~ "Microcoleus", 
                                     "non_e_diatoms" ~ "Non-Epithemia Diatoms", 
                                     "nostoc" ~ "Nostoc", 
                                     "oscillatoria" ~ "Oscillatoria",
                                     "other_coccoids" ~ "Other Coccoids", 
                                     "rare" ~ "Rare"))


#Standardize TAC mats to Anabaena abundances
microscopy_TAC_stand <- microscopy %>% 
  dplyr::filter(sample_type == "TAC") %>% 
  mutate(across(anabaena_and_cylindrospermum:rare, ~ . / anabaena_and_cylindrospermum)) %>% 
  select(!anabaena_and_cylindrospermum) %>% 
  mutate(total = rowSums(select(., e_diatoms:rare))) %>% 
  mutate(across(e_diatoms:rare, ~. / total)) %>% 
  select(!total) %>% 
  pivot_longer(e_diatoms:rare, names_to = "Species", values_to = "Abundance")


#Standardize TM mats to Microcoleus abundances
microscopy_TM_stand <- microscopy %>% 
  dplyr::filter(sample_type == "TM") %>% 
  mutate(across(anabaena_and_cylindrospermum:rare, ~ . / microcoleus)) %>% 
  select(!microcoleus) %>% 
  mutate(total = rowSums(select(., anabaena_and_cylindrospermum:rare))) %>% 
  mutate(across(anabaena_and_cylindrospermum:rare, ~. / total)) %>% 
  select(!total) %>% 
  pivot_longer(anabaena_and_cylindrospermum:rare, names_to = "Species", values_to = "Abundance")



#------Plots and Graphs---------------------------------------

#Create a color palette
mycols <- brewer.pal(n = 11, name = "Set3")
mypal <- palette(mycols)
names(mypal) = c("Anabaena", "Epithemia Diatoms", "Geitlerinema", 
                 "Green Algae", "Leptolyngbya", "Microcoleus", 
                 "Non-Epithemia Diatoms", "Nostoc", "Oscillatoria",
                 "Other Coccoids", "Rare")
colScale <- scale_fill_manual(name = "Species", values = mypal)

#Plot timeseries of within-mat species per year
ggplot(subset(microscopy, sample_type %in% "TM"), aes(x = field_date, y = Abundance, fill = Species)) +
  facet_grid(reach ~ year, scales = "free") +
  geom_col(position = "fill") +
  theme(legend.position="bottom") + 
  labs(title = "Target Microcoleus", x = "Date", y = "Proportion of Relative Abundance") +
  colScale

#Plot timeseries, but remove Microcoleus
ggplot(subset(microscopy, sample_type %in% "TM" & !(Species %in% "Microcoleus")), 
       aes(x = field_date, y = Abundance, fill = Species)) +
  facet_grid(reach ~ year, scales = "free") +
  geom_col(position = "fill") + 
  theme(legend.position="bottom") + 
  labs(title = "Target Microcoleus", x = "Date", y = "Proportion of Relative Abundance") +
  colScale

#Plot timeseries of within-mat species per year
ggplot(subset(microscopy, sample_type %in% "TAC"), aes(x = field_date, y = Abundance, fill = Species)) +
  facet_grid(reach ~ year, scales = "free") +
  geom_col(position = "fill") +
  theme(legend.position="bottom") + 
  labs(title = "Target Anabaena", x = "Date", y = "Proportion of Relative Abundance") +
  colScale

#Plot timeseries, but remove Anabaena
ggplot(subset(microscopy, sample_type %in% "TAC" & !(Species %in% "Anabaena")), 
       aes(x = field_date, y = Abundance, fill = Species)) +
  facet_grid(reach ~ year, scales = "free") +
  geom_col(position = "fill") + 
  theme(legend.position="bottom") + 
  labs(title = "Target Anabaena", x = "Date", y = "Proportion of Relative Abundance") +
  colScale

