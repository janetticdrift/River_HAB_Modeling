#Code to model within-mat dynamics
library(tidyverse)

withinmat <- read.csv(here::here("data/16S_2024.csv")) #2024 data

#Must remove the percentages in percent_comp
withinmat <- withinmat %>% 
  mutate(across(everything(), str_remove_all, "%"))
withinmat$percent_comp <- as.numeric(as.character(withinmat$percent_comp))
