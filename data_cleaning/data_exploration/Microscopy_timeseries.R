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
  filter(grepl("SFE", site_reach)) %>%  #Keep sites that include string "SFE"
  separate(site_reach, into=c("site", "location", "reach"), 
           sep="-") %>% #Split location columns into separate categories
  filter(grepl("M", location)) %>% #Keep Miranda sites, remove Standish-Hicky (SH) sites
  mutate(field_date = as.Date(field_date, format="%m/%d/%y")) %>%  #Fix date data type
  mutate(year = year(field_date)) %>%  #create year column
  relocate(year, .after = field_date) %>%   #reorganize column order
  select(!non_algal) %>% #Remove column measuring sediment amount
  pivot_longer(microcoleus:aphanothece, names_to = "Species", values_to = "Abundance")    #Change layout of dataframe
  
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
  #Bind together dataframes
  microscopy <- rbind(averaged_slides, processed_slides)
  
  
#Do TAC samples all actually contain Anabaena?
TAC <- microscopy %>% 
  filter(sample_type == "TAC" & Species == "anabaena_and_cylindrospermum")

which(TAC$Abundance < 10) #Are there any instances where Abundance is small? FALSE = no
which(TAC$Abundance < 20)
  #Some instances where Anabaena is less than 20% in TAC sample

#Are there any species that never appear?
non_occurences <- microscopy %>% 
  group_by(year, sample_type, Species) %>% 
  dplyr::summarise(Sum = sum(Abundance))

rownum <- which(non_occurences$Sum <= 0) #Which samples only find 1% summed across all sampled dates per year

rare_species <- non_occurences %>% 
  #dplyr::filter(year == "2022") %>% 
  ungroup() %>% 
  dplyr::slice(rownum)

rare_names <- unique(rare_species$Species) #Names of the species that are rare

#Remove non-occuring species, and combine those that are super rare
microscopy1 <- microscopy %>% 
  #dplyr::filter(!Species == "tolypothrix") %>% #Remove non-occuring
  dplyr::mutate(Species = fct_collapse(Species,
                                        rare = c(rare_names)))


#Standardize TAC mats to Anabaena abundances
TAC_ref_group <- microscopy1 %>% 
  dplyr::filter(sample_type == "TAC") %>% 
  dplyr::filter(Species == "anabaena_and_cylindrospermum") %>%
  group_by(year) %>% 
  dplyr::summarise(ref_mean = mean(Abundance),
            ref_sd = sd(Abundance))  #Calculate mean and SD for Anabaena per year

#Standardize all TAC mat species using reference group
microscopy_TAC_stand <- microscopy1 %>% 
  dplyr::filter(sample_type == "TAC") %>% 
  mutate(stand_abund = (Abundance - TAC_ref_group$ref_mean) / TAC_ref_group$ref_sd) 



#------Plots and Graphs---------------------------------------




#Plot timeseries of within-mat species per year
ggplot(subset(microscopy1, sample_type %in% "TM"), aes(x = field_date, y = Abundance, fill = Species)) +
  facet_wrap(~year, scales = "free") +
  geom_col(position = "fill") + #Can't see the detail well, subset per year
  theme(legend.position="bottom") + 
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Target Microcoleus")

ggplot(subset(microscopy1, sample_type %in% "TAC"), aes(x = field_date, y = Abundance, fill = Species)) +
  facet_wrap(~year, scales = "free") +
  geom_col(position = "fill") + #Can't see the detail well, subset per year
  theme(legend.position="bottom")+
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Target Anabaena")

#Coding notes: add year %in% 2022 to subset if splitting up years again

#2022
#Histogram of Microcolus in TM samples
ggplot(subset(microscopy1, year %in% 2022 & 
                Species %in% "microcoleus" & 
                sample_type %in% "TM"), 
       aes(x = Abundance, fill = Species)) +
  geom_histogram(binwidth = 1, position="dodge") +
  labs(title = "2022 Counts") +
  theme(legend.position="bottom")

#Histogram of Anabaena in TA samples
ggplot(subset(microscopy, year %in% 2022 & 
                Species %in% "anabaena_and_cylindrospermum" &
                sample_type %in% "TAC"), 
       aes(x = Abundance, fill = Species)) +
  geom_histogram(binwidth = 1, position="dodge") +
  labs(title = "2022 Counts") +
  scale_fill_manual(values = c("anabaena_and_cylindrospermum" = "#00BFC4")) +
  theme(legend.position="bottom")



#2023

#Histogram of Microcolus in TM samples
  ggplot(subset(microscopy, year %in% 2023 & 
                  Species %in% "microcoleus" & 
                  sample_type %in% "TM"), 
         aes(x = Abundance, fill = Species)) +
    geom_histogram(binwidth = 1, position="dodge") +
    labs(title = "2023 Counts") +
    theme(legend.position="bottom")
  
  #Histogram of Anabaena in TA samples
  ggplot(subset(microscopy, year %in% 2023 & 
                  Species %in% "anabaena_and_cylindrospermum" &
                  sample_type %in% "TAC"), 
         aes(x = Abundance, fill = Species)) +
    geom_histogram(binwidth = 1, position="dodge") +
    labs(title = "2023 Counts") +
    scale_fill_manual(values = c("anabaena_and_cylindrospermum" = "#00BFC4")) +
    theme(legend.position="bottom")
  

  
  
#Plots using standardized values
  ggplot(subset(microscopy_TAC_stand), aes(x = field_date, y = stand_abund, color = Species)) +
    facet_wrap(~year, scales = "free") +
    geom_point() +
    theme(legend.position="bottom")+
    scale_fill_brewer(palette = "Set3") +
    labs(title = "Target Anabaena")
