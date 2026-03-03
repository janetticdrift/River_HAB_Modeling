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

pcr_atx <- read.csv(here::here("data/qPCR_Anatoxins.csv")) %>% 
  dplyr::filter(!grepl("analysis", UNR_Notes)) %>%  #Remove samples noted as not for analysis ("floating" and "gravel")
  rename(field_date = date) %>% 
  dplyr::mutate(field_date = as.Date(field_date, format="%m/%d/%y")) %>%  # Fix date data type
  dplyr::mutate(
    across(c(ana_C.uL:dhATXa_ug_g),
           ~ case_when(. == "bdl" ~ '0',
                       . != "bdl" ~ .))) %>% # Convert non-detects to zeros
  dplyr::mutate(
    across(c(ana_C.uL, ana_C.uL_rerun,
             nif.uL, nif.uL_rerun),
           ~ log10(.x + pseudocount))  #Add small pseudocount to all copy numbers, then log transform
  ) %>% 
  pivot_longer(   #Pivot dataframe to create a group for run vs. rerun in the next function
    cols = c(ana_C.uL, ana_C.uL_rerun,
             nif.uL, nif.uL_rerun),
    names_to = "gene_run",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    gene = case_when(                     #Remove extra characters from gene name category names
      str_detect(gene_run, "ana_C") ~ "anaC",
      str_detect(gene_run, "nif")  ~ "nif"
    ),
    run_type = if_else(str_detect(gene_run, "rerun"), "rerun", "original") #Create run_type grouping
  ) %>%
  select(-gene_run) %>% #Remove redundant column
  pivot_wider(names_from = gene, values_from = value) %>%
  dplyr::mutate(
    across(c(ATX_all_ug_g, ATXa_ug_g,    #Convert ATX values from character to numeric
             HTXa_ug_g, dhATXa_ug_g),
           ~ as.numeric(.))
  )



#------------------------------------------------------------------------------
################
#GRAPHING GENES#
################

#nif ~ anaC, exclude zeros from regression line
ggplot(pcr_atx, aes(x = anaC, y = nif, color = mat)) +
  facet_wrap(~year, scales = "free") +  #Split plot between the two years, scales = "free" means each year has its own Y axis 
  geom_point(aes(shape = run_type), size = 2) +
  geom_smooth(data = . %>% filter(anaC != -4 & nif != -4), # filter out pseudo "0" values
              method = "lm", se = FALSE) +
  stat_regline_equation(                  #Calculate r^2 values
    data = . %>% filter(anaC != -4 & nif != -4),   #Remove pseudo "0" values from calculating r^2
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    parse = TRUE,
    label.x.npc = "left",
    label.y.npc = "top") +
  labs(x = "Log anaC gene (copies/ng)", y = "Log nif gene (copies/ng)", title = "anaC and nif") +
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

