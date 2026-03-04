############################################
###Anatoxin Data Cleaning and Exploration###
############################################

#Packages----
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#Read in datasheet and clean data

pseudocount <- 1e-4 #Set a small number for 0s, so as to not return error when logging

pcr_atx_combine <- read.csv(here::here("data/qPCR_Anatoxins.csv")) %>% 
  pivot_longer(cols = c(sample:nif.uL_rerun), 
               names_to = c("measure", ".value"), 
               names_pattern = "(.*)_(.*)") %>%
  # Combine the rows where duplicates exist based on unique ID (field_date, year, site, mat)
  group_by(field_date, year, site, mat) %>%
  # Create a new row for duplicate samples, add "qPCR duplicate" to Notes if 2 samples exist
  mutate(Notes = if_else(n() > 1, "qPCR duplicate", NA_character_)) %>%
  ungroup() %>%
  # Sort by the measure and make sure NAs are handled as necessary (e.g., replace NA with rerun values if needed)
  arrange(field_date, year, site, mat)

  
  
  dplyr::filter(!grepl("analysis", Notes)) %>%  #Remove samples not for analysis ("floating" and "gravel")
  rename(field_date = date) %>% 
  dplyr::mutate(field_date = as.Date(field_date, format="%m/%d/%y")) %>%  #Fix date data type
  dplyr::mutate(ana_C.uL = ifelse(UNR_Notes == "PCR duplicate", NA, ana_C.uL)) %>% #Remove 2024 first analyses if they were questionable duplicates
  dplyr::mutate(
    # # Normalized anaC gene copies (copies per ng DNA)
    # normalized_anaC = ana_C.uL / DNA_conc_ng_uL,
    # normalized_anaC_rerun = ana_C.uL_rerun / DNA_conc_rerun,
    # # Normalized nif gene copies (copies per ng DNA)
    # normalized_nif = nif.uL / DNA_conc_ng_uL,
    # normalized_nif_rerun = nif.uL_rerun / DNA_conc_rerun,
    
    # Set any "bdl" values to 0
    across(c(ana_C.uL:dhATXa_ug_g),
           ~ case_when(. == "bdl" ~ '0',
                       . != "bdl" ~ .)),
    
    # Convert values from character to numeric
    across(c(ana_C.uL:dhATXa_ug_g), ~ as.numeric(.)),
    
    # Log10 gene counts
    across(c(ana_C.uL, ana_C.uL_rerun,
             nif.uL, nif.uL_rerun),
           ~ log10(.x + pseudocount)),
    
    # Calculate ATX normalized with Ash-Free Dry Mass
    norm_ATX_all_ug_g = ATX_all_ug_g / org_matter_percent
    
  ) %>% 
  pivot_longer(   #Pivot dataframe to create a group for run vs. rerun in the next function
    cols = c(ana_C.uL, ana_C.uL_rerun,
             nif.uL, nif.uL_rerun),
    names_to = "gene_run",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    gene = case_when(           #Remove extra characters from gene name category names
      str_detect(gene_run, "ana_C") ~ "anaC",
      str_detect(gene_run, "nif")  ~ "nif"),
    run_type = if_else(str_detect(gene_run, "rerun"), "rerun", "original")) %>% #Create run_type grouping
  dplyr::select(-gene_run) %>% 
  pivot_wider(names_from = gene, values_from = value)


#This dataframe has non-normalized gene copy values, and normalized ATX values


#------------------------------------------------------------------------------
################
#GRAPHING GENES#
################

#nif ~ anaC, exclude zeros from regression line
ggplot(subset(pcr_atx, year %in% 2024), aes(x = anaC, y = nif, color = mat)) +
  facet_wrap(~year, scales = "free") +  #Split plot between the two years, scales = "free" means each year has its own Y axis 
  geom_point(aes(shape = run_type), size = 2) 
  geom_smooth(data = . %>% dplyr::filter(anaC != -4 & nif != -4), # filter out pseudo "0" values
              method = "lm", se = FALSE) +
  stat_regline_equation(                  #Calculate r^2 values
    data = . %>% dplyr::filter(anaC != -4 & nif != -4),   #Remove pseudo "0" values from calculating r^2
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    label.x.npc = "left",
    label.y.npc = "top",
    show.legend = FALSE) +
  labs(x = "Log anaC gene (copies/ng)", y = "Log nif gene (copies/ng)") +
  theme_bw()

#nif ~ anaC, include zeros in regression line
ggplot(pcr_atx, aes(x = anaC, y = nif, color = mat)) +
  facet_wrap(~year, scales = "free") +
  geom_point(aes(shape = run_type), size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_regline_equation(
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    label.x.npc = "left",
    label.y.npc = "top") +
  labs(x = "anaC gene (copies/ng)", y = "nif gene (copies/ng)") +
  theme_bw()

##############
#GRAPHING ATX#
##############

#ana_C vs Anatoxins, exclude zeros from regression line
ggplot(pcr_atx, aes(x = anaC, y = log(ATX_all_ug_g), color = mat)) + #ATXs log-transformed
  facet_wrap(~year, scales = "free") +
  geom_point( size = 2) +
  geom_smooth(data = . %>% filter(anaC != -4 & ATX_all_ug_g != 0),
              method = "lm", se = FALSE) +
  stat_regline_equation(
    data = . %>% filter(anaC != -4 & ATX_all_ug_g != 0),
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    label.x.npc = "left",
    label.y.npc = "top") +
  labs(x = "Log anaC gene (copies/ng)", y = "Log total anatoxin (ug/g)", title = "ATX and anaC") +
  scale_colour_brewer(type = "qual") +    #Use different colors for anatoxin plot
  theme_bw()

