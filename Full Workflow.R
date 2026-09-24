###########################
#Analysis Methods Workflow
###########################
#This file consolidates all modeling and analytical workflows and sources them directly 
#here. This will allow the user to run each step of the paper's analtyical methodology step
#by step, from estimating algal abundance latent states, to predicting algal abundances,
#to estimating toxin latent states using algal predictions, to predicting toxins. Also 
#included is the analysis on algal sensitivity to different environmental parameters,
#and calculations of model fit metrics and evaluation of predictive performance.

#The primary purpose of this file is to streamline the methodology's workflow from several
#inter-dependent files into one clean process. It also serves as a Table of Contents,
#listing the steps taken to conduct the analyses, and the filed used and their locations in
#each folder.

#Sourcing entire files means that the files containing Bayesian models may take several 
#minutes to run first before figures can be generated.

library(here)
library(patchwork)

######Acronym Key######
#RW = River-Wide
#WM = Within-Mat
#TM = Target Microcoleus
#TAC = Target Anabaena/Cylindrospermum 


###################################################################
###Benthic Percent Cover and Microscopy Bayesian Models
###################################################################

#Runs the Bayesian models
source(here::here("data_analysis/Running Stan Models/Latent_States_Models.R"))

  #This model also contains code for running ELPD (loo package) and DIC. These calculations
#remain in this file as they requires the full model fit, which are too large to save
#as RDS files.

m1.1.elpd
m1.2.elpd
m1.3.elpd
m1.4.elpd
m1.5.elpd
m2.elpd

#Compares observed data to latent states
#Percent cover model
source(here::here("data_analysis/Compare Obs Vs Modeled Outputs/River_Wide_model_vs_real.R"))
#Microscopy model
source(here::here("data_analysis/Compare Obs Vs Modeled Outputs/Within_Mat_model_vs_real.R"))

#----------Relevant figures----------#

#Supplemental
obs.v.real_RW_all
obs.v.real_RW_biotic
obs.v.real_RW_abiotic
obs.v.real_RW_abioticnonut
obs.v.real_RW_trueabiotic

#Specify the locations of plots in the layout. 2 columns, 6 rows. J appears twice to take
  #up two row spaces in the bottom of the right column.
design <- "
AB
CD
EF
GH
IJ
KJ
"

#Wrap_plots lists plots in the order of left->right, top->bottom
SupFigure1 <- wrap_plots(p1,p3,p2,p4,p5,p7,p6,p8,p9,
  guide_area(),   # collected legend
  p10,
  design = design) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "horizontal")

pdf(file = "Figures/SupFigure1.pdf", width = 10, height = 6)
(SupFigure1)
dev.off()

obs.v.real_WM_all

pdf(file = "Figures/SupFigure2.pdf", width = 10, height = 5)
(obs.v.real_WM_all)
dev.off()

###################################################################
###Benthic Percent Cover and Microscopy Prediction Figures
###################################################################

#Percent-cover predictions
source(here::here("data_analysis/Running Predictions/Predictions_River_Wide.R"))
#Microscopy predictions
source(here::here("data_analysis/Running Predictions/Predictions_Within_Mat.R"))

#----------Relevant figures----------#

Figure2 <- (p2.breaks / p4.breaks / p8.breaks) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "vertical") 

pdf(file = "Figures/Figure2.pdf", width = 10, height = 5)
(Figure2)
dev.off()

#Supplemental 

SupFigure3 <- wrap_plots(p1,p3,p2.complete,p4.complete,p5,p7,p6.complete,p8.complete,p9,
                         guide_area(),   # collected legend
                         p10.complete,
                         design = design) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "right", legend.box = "horizontal")

pdf(file = "Figures/SupFigure3.pdf", width = 10, height = 5)
(SupFigure3)
dev.off()

RWplot.all
RWplot.biotic
RWplot.abiotic
RWplot.abioticnonut
RWplot.trueabiotic

WMplot.all

###################################################################
###Sensitivity Analyses: Equilibrium Abundances
###################################################################

source(here::here("data_analysis/Sensitivity Analysis.R"))

Figure3 <- (micro_eq /ana_eq) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect", axis_titles = "collect_x")

pdf(file = "Figures/Figure3.pdf", width = 10, height = 5)
(Figure3)
dev.off()

#Supplemental 

pdf(file = "Figures/SupFigure4.pdf", width = 10, height = 5)
(micro_total_eq)
dev.off()

pdf(file = "Figures/SupFigure5.pdf", width = 10, height = 5)
(ana_total_eq)
dev.off()


###################################################################
###Toxin Concentration Bayesian Models
###################################################################

#Runs the Bayesian models
source(here::here("data_analysis/Running Stan Models/Toxin Models.R"))
#Compares observed data to latent states
source(here::here("data_analysis/Compare Obs Vs Modeled Outputs/Toxins_model_vs_real.R"))

#----------Relevant figures----------#

#Supplemental

pdf(file = "Figures/SupFigure6.pdf", width = 10, height = 5)
(lagplot)
dev.off()


design <- "
AB
CD
EF
"

SupFigure7 <- wrap_plots(obs.v.real.plots$TOX_RW_TM_All,obs.v.real.plots$TOX_RW_TM_Biotic,
                         obs.v.real.plots$TOX_RW_TM_Abiotic,obs.v.real.plots$TOX_RW_TM_AbioticNoNut,
                         obs.v.real.plots$TOX_RW_TM_TrueAbiotic,
                         guide_area(),   # collected legend
                         design = design) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axes = "collect_y") &
  theme(legend.position = "right", legend.box = "horizontal")

pdf(file = "Figures/SupFigure7.pdf", width = 10, height = 6)
(SupFigure7)
dev.off()

SupFigure8 <- wrap_plots(obs.v.real.nolag.plots$TOX_RW_TM_All,obs.v.real.nolag.plots$TOX_RW_TM_Biotic,
                         obs.v.real.nolag.plots$TOX_RW_TM_Abiotic,obs.v.real.nolag.plots$TOX_RW_TM_AbioticNoNut,
                         obs.v.real.nolag.plots$TOX_RW_TM_TrueAbiotic,
                         guide_area(),   # collected legend
                         design = design) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axes = "collect_y") &
  theme(legend.position = "right", legend.box = "horizontal")

pdf(file = "Figures/SupFigure.pdf", width = 10, height = 6)
(SupFigure8)
dev.off()


obs.v.real_TOX_RW_TAC

pdf(file = "Figures/SupFigure9.pdf", width = 10, height = 6)
(obs.v.real_TOX_RW_TAC)
dev.off()

obs.v.real_TOX_WM_TM

pdf(file = "Figures/SupFigure10.pdf", width = 10, height = 6)
(obs.v.real_TOX_WM_TM)
dev.off()

###################################################################
###Toxin Concentration Prediction Figures
###################################################################

source(here::here("data_analysis/Running Predictions/Predictions_Toxins.R"))

#----------Relevant figures----------#

Figure4 <- (all_model_plots$All / all_model_plots$Biotic / all_model_plots$AbioticNoNut / envplot) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")

#OR#

Figure4 <- (all_model_plots$All / envplot / nolag_model_plots$All) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")

pdf(file = "Figures/Figure4.pdf", width = 10, height = 5)
(Figure4)
dev.off()

#Supplement
SupFigure10 <- wrap_plots(all_model_plots$All,all_model_plots$Biotic,
                         all_model_plots$Abiotic,all_model_plots$AbioticNoNut,
                         all_model_plots$TrueAbiotic,
                         guide_area(),   # collected legend
                         design = design) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axes = "collect_y") &
  theme(legend.position = "right", legend.box = "horizontal")

SupFigure11 <- wrap_plots(nolag_model_plots$All,nolag_model_plots$Biotic,
                          nolag_model_plots$Abiotic,nolag_model_plots$AbioticNoNut,
                          nolag_model_plots$TrueAbiotic,
                          guide_area(),   # collected legend
                          design = design) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axes = "collect_y") &
  theme(legend.position = "right", legend.box = "horizontal")

RWToxplot.TAC
WMToxplot


###################################################################
###Model Fit Indices
###################################################################

source(here::here("data_analysis/Model Fit Calculations.R"))

#----------Relevant figures----------#

Figure5 <- (percover.metricplot | tox.metricplot) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal")


