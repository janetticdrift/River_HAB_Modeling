########
##Learning to use gjamTIME with SFE Data##

#This code works through the gjamTime function, an extension of the gjam package, using
#algae bloom and toxin data from the South Fork Eel River in California, USA.

#The experimental setup includes...


########

#Library
library(plyr)
library(tidyverse)
library(ggpubr)
library(gjam)
library(devtools)
library(here)

#Source the gjamTIME functions
d <- "https://github.com/jimclarkatduke/gjam/blob/master/gjamTimeFunctions.R?raw=True"
source_url(d)

#Read in data
algae <- read.csv(here::here("data/percover_byreach.csv"))
ATX <- read.csv(here::here("data/cyano_atx.csv"))


#Needed dataframes
#####################
#xdata

xdata.temp <- algae %>% 
  dplyr::filter(site %in% "SFE-M") %>% 
  dplyr::filter(grepl('2022', field_date)) %>% 
  dplyr::select(!c(site, reach, micro_transects:proportion_ana_cyl_transects)) %>% 
  tidyr::pivot_longer(cols = green_algae:other_nfixers, names_to = "Species", 
                      values_to = "Abundance") %>% 
  dplyr::mutate(intercept = 1)

xdata <- xdata.temp %>% 
  dplyr::filter(Species == "green_algae") %>% #Remove other species, so there is only one plot entry per year
  dplyr::select(c(field_date, site_reach, intercept)) %>% 
  dplyr::mutate(field_date = as.numeric(1:length(field_date))) %>% 
  as.data.frame()
  #dplyr::left_join(gjamclimdat, by = "field_date?") 
  
#####################
#ydata

ydata <- xdata.temp %>% 
  tidyr::pivot_wider(names_from = Species, values_from = Abundance) %>% #Pivot data to put species along the top and populate matrix with abundances
  dplyr::select(!(field_date:intercept)) %>% #"Deselect" year, plot, and intercept column from the dataframe
  dplyr::mutate_all(funs(./100, )) %>% 
  as.matrix() #Turn it into a matrix

#####################
#edata

edata <- ydata
edata[edata>-1] <- 1  

#####################
#Formula (1 parameter(s))

formula <- ~ 1

#####################
#Fill in missing values

timeCol   <- 'field_date'
groupCol  <- 'site_reach'
groupVars <- c( 'site_reach' ) #identifies those columns in xdata that are fixed for the group, so they can be filled in
fillmissdata <- gjamFillMissingTimes(xdata, ydata, edata, groupCol, timeCol,
                                     FILLMEANS = T, groupVars = groupVars,
                                     typeNames = 'FC', missingEffort = .1)   #FC = Fractional composition

xdatamiss  <- fillmissdata$xdata
ydatamiss  <- fillmissdata$ydata
edatamiss  <- fillmissdata$edata
tlist  <- fillmissdata$timeList #bookkeeping objects used for vectorized operations
snames <- colnames(ydata) #species names
effort <- list(columns = 1:length(unique(xdata.temp$Species)), values = edatamiss) #Fill in observation effort

#####################
#Set priors

#rhoPrior is a list indicating lo and hi observed abundance change values (changes per time increment).
#plus or minus 50% growth rate per time increment
rhoPrior  <- list(lo = list(intercept = -.5), 
                  hi = list(intercept = .5))


#Make a priorList for gjamTimePrior
priorList <- list(formulaRho = as.formula(~ 1),  
                      rhoPrior = rhoPrior)

#Organize prior parameter values into a list
tempdat <- gjamTimePrior(xdatamiss, ydatamiss, edatamiss, priorList)
timeList <- mergeList(tlist, tempdat) #Update the list

#####################
#Model fitting

modelList <- list(typeNames = 'FC', ng = 4000, burnin = 1000,  #FC = Fractional composition
                      timeList = timeList, effort = effort) 

outputAR <- gjam(formula, xdata=xdatamiss, ydata=ydatamiss, modelList=modelList)

#Plot
plotPars  <- list(PLOTALLY=T, 
                  SAVEPLOTS = F)

gjamPlot(outputAR, plotPars)
