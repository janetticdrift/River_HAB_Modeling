#Code to model within-mat dynamics
library(tidyverse)

withinmat <- read.csv(here::here("data/16S_2024.csv")) #2024 data
genecopy <- read.csv(here::here("data/qPCR_2024.csv")) #2024 data

#Clean data: Must remove the percentages in percent_comp
withinmat <- withinmat %>% 
  mutate(across(everything(), str_remove_all, "%"))
withinmat$percent_comp <- as.numeric(as.character(withinmat$percent_comp))

library(plotly)
library(ggplot2)


p1 <- ggplot(genecopy, aes(x = DNA_conc, y = ana_C.uL, color = mat, label = sample)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 2000)) +
  labs(title = "anaC Gene Copies", x = "DNA concentration ng/uL", 
       y = "anaC Copy Count")

p2 <- ggplot(genecopy, aes(x = DNA_conc, y = nif.uL, color = mat, label = sample)) +
  geom_point() +
  coord_cartesian(ylim = c(0, 2000)) +
  labs(title = "nif Gene Copies", x = "DNA concentration ng/uL", 
       y = "nif Copy Count")

ggarrange(
  p1, p2, labels = c("A", "B"),
  common.legend = TRUE, legend = "bottom"
)


ggplotly(plotgenecopy, tooltip = "sample")
