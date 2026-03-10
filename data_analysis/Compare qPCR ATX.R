############################################
###Anatoxin Data Cleaning and Exploration###
############################################

#Packages----
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(here)

#Read in datasheet and clean data

pseudocount <- 1e-2 #Set a small number for 0s, so as to not return error when logging

pcr_atx_df <- read.csv(here::here("data/qPCR_Anatoxins.csv")) %>% 
  dplyr::filter(!grepl("analysis", UNR_Notes)) %>%  #Remove samples not for analysis ("floating" and "gravel")
  dplyr::rename(field_date = date) %>% 
  dplyr::mutate(field_date = as.Date(field_date, format="%m/%d/%y")) %>%  #Fix data type
  mutate(samp_rerun = as.character(samp_rerun),
         sample = ifelse(sample == "", NA, sample),
         ana_C.uL = ifelse(sample == "", NA, ana_C.uL),
         nif.uL = ifelse(sample == "", NA, nif.uL)). #Fill in NAs for samples that were not collected


pcr_combine <- pcr_atx_df %>% 
  # Create two datasets: "original" and "rerun"
  dplyr::select(field_date, year, site, mat,
                sample, DNA_conc_ng_uL, ana_C.uL, nif.uL, 
                ATX_all_ug_g:org_matter_percent, UNR_Notes) %>%
  dplyr::mutate(run_type = "original") %>%
  bind_rows(
    pcr_atx_df %>%
      dplyr::select(field_date, year, site, mat,
                    samp_rerun, DNA_conc_rerun, ana_C.uL_rerun, nif.uL_rerun, 
                    ATX_all_ug_g:org_matter_percent, UNR_Notes) %>% 
      dplyr::mutate(run_type = "rerun",
             sample = samp_rerun,
             DNA_conc_ng_uL = DNA_conc_rerun,
             ana_C.uL = ana_C.uL_rerun,
             nif.uL = nif.uL_rerun)) %>%
  dplyr::filter(!(is.na(sample) & is.na(DNA_conc_ng_uL) & is.na(ana_C.uL) & is.na(nif.uL))) %>% #Find duplicates by removing NA rows
  group_by(field_date, year, site, mat) %>%
  dplyr::mutate(UNR_Notes = ifelse(n() > 1, "qPCR duplicate", NA)) %>% #Note duplicates (grouped by the above cols)
  ungroup() %>% 
  dplyr::select(1:15) #Keep combined columns

#Vector for filtering out sus samples
sus_anaC <- c("12", "13", "32", "33", "34", "35")
sus_nif <- c("45", "39")
  
pcr_atx <- pcr_combine %>% 
  dplyr::mutate(
    # # Normalized anaC gene copies (copies per ng DNA)
    # normalized_anaC = ana_C.uL / DNA_conc_ng_uL,
    # # Normalized nif gene copies (copies per ng DNA)
    # normalized_nif = nif.uL / DNA_conc_ng_uL,
    
    #Note that some of the rerun samples did not have an updated DNA_conc value
    
    # Set any "bdl" values to 0
    across(c(ana_C.uL:dhATXa_ug_g),
           ~ case_when(. == "bdl" ~ '0',
                       . == "nd" ~ '0',
                       . != "bdl" ~ .,
                       . != "nd" ~ .)),
    
    # Convert values from character to numeric
    across(c(ana_C.uL:dhATXa_ug_g), ~ as.numeric(.)),
    
    # Log10 gene counts
    across(c(ana_C.uL, nif.uL),
           ~ log10(.x + pseudocount)),
    
    # Calculate ATX normalized with Ash-Free Dry Mass
    norm_ATX_all_ug_g = (ATX_all_ug_g / org_matter_percent) + pseudocount
  )  %>% 
  dplyr::filter(!c(year == "2024" & sample %in% sus_anaC),
         !c(year == "2024" & sample %in% sus_nif))


#This dataframe has non-normalized gene copy values, and normalized ATX values


#------------------------------------------------------------------------------
################
#GRAPHING GENES#
################

#nif ~ anaC, exclude zeros from regression line
ggplot(pcr_atx, aes(y = ana_C.uL, x = nif.uL, color = mat)) +
  facet_wrap(~year, scales = "free") +  #Split plot between the two years, scales = "free" means each year has its own Y axis 
  geom_point(aes(shape = run_type), size = 2) +
  geom_smooth(data = . %>% dplyr::filter(ana_C.uL != -2 & nif.uL != -2), # filter out pseudo "0" values
              method = "lm", se = FALSE) +
  stat_regline_equation(                  #Calculate r^2 values
    data = . %>% dplyr::filter(ana_C.uL != -2 & nif.uL != -2),   #Remove pseudo "0" values from calculating r^2
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    label.x.npc = "left",
    label.y.npc = "top",
    show.legend = FALSE) +
  labs(y = "Log10 anaC gene (copies/ng)", x = "Log10 nif gene (copies/ng)") +
  theme_bw()

#nif ~ anaC, include zeros in regression line
ggplot(pcr_atx, aes(x = ana_C.uL, y = nif.uL, color = mat)) +
  facet_wrap(~year, scales = "free") +
  geom_point(aes(shape = run_type), size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_regline_equation(
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    label.x.npc = "left",
    label.y.npc = "top",
    show.legend = F) +
  labs(x = "anaC gene (copies/ng)", y = "nif gene (copies/ng)") +
  theme_bw()

##############
#GRAPHING ATX#
##############

#ana_C vs Anatoxins, exclude zeros from regression line
ggplot(pcr_atx, aes(x = ana_C.uL, y = log(ATX_all_ug_g), color = mat)) + #ATXs log-transformed
  facet_wrap(~year, scales = "free") +
  geom_point( size = 2) +
  geom_smooth(data = . %>% filter(ana_C.uL != -2 & ATX_all_ug_g != 0),
              method = "lm", se = FALSE) +
  stat_regline_equation(
    data = . %>% filter(ana_C.uL != -2 & ATX_all_ug_g != 0),
    aes(label = paste(..rr.label..)),
    formula = y ~ x,
    label.x.npc = "left",
    label.y.npc = "top",
    show.legend = F) +
  labs(x = "Log anaC gene (copies/ng)", y = "Total anatoxin (ug/g)", title = "ATX and anaC") +
  scale_colour_brewer(type = "qual") +    #Use different colors for anatoxin plot
  theme_bw()

