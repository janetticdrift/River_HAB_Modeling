#Code to model within-mat dynamics
library(tidyverse)

#Microscopy
withinmat <- read.csv(here::here("data/16S_2024.csv")) #2024 data
#Gene copies
genecopy <- read.csv(here::here("data/qPCR_2024.csv")) #2024 data
#Ash Free Dry Mass
afdm <- read.csv(here::here("data/Ashfree Dry Mass 2024.csv"))

#Cleaning data: Remove the percentages in percent_comp
withinmat <- withinmat %>% 
  mutate(across(everything(), str_remove_all, "%"))
withinmat$percent_comp <- as.numeric(as.character(withinmat$percent_comp))

#Cleaning data: Organize AFDM, remove improbable row
afdm <- afdm %>% 
  dplyr::filter(River %in% 'South Fork Eel')

library(plotly)
library(ggplot2)


p1 <- ggplot(genecopy, aes(x = DNA_conc, y = log(ana_C.uL+1), color = mat, label = sample)) +
  geom_point() +
  #coord_cartesian(ylim = c(-50, 20)) +
  labs(title = "anaC Gene Copies", x = "DNA concentration ng/uL", 
       y = "anaC Copy Count")

p2 <- ggplot(genecopy, aes(x = DNA_conc, y = log(nif.uL+1), color = mat, label = sample)) +
  geom_point() +
  #coord_cartesian(ylim = c(0, 2000)) +
  labs(title = "nif Gene Copies", x = "DNA concentration ng/uL", 
       y = "nif Copy Count")

p3 <- ggplot(genecopy, aes(x = log(ana_C.uL+1), y = log(nif.uL+1), color = mat, label = sample)) +
  geom_point() +
  #coord_cartesian(ylim = c(0, 2000)) +
  labs(title = " Gene Copies", x = "ana_C copy count", 
       y = "nif Copy Count")

ggarrange(
  p1, p2, labels = c("A", "B"),
  common.legend = TRUE, legend = "bottom"
)


#anaC gene copy count by DNA concentration, log-transformed
ggplotly(p1, tooltip = "sample")
#nif gene copy count by DNA concentration, log-transformed
ggplotly(p2, tooltip = "sample")
#gene copy count in Ana vs Micro
ggplotly(p3, tooltip = "sample")


