#####
#HAB Model
#####

#Set up a for loop that runs through each species (1 to start) in the model. This model
#uses theoretical numbers, no real data
#Set up environmental effect on DD DI  - how will that later translate to real env effects?
  #See Kleinhesselink 2018

#Library
library(tidyverse)

########
#Single species
########
time <- 50

  #Set state variables
  N1 <- rep(NA,time)
  #Starting abundance
  N1[1] <- 2#log(50) 

  #Set parameters
  r <- 1 #runif(1, .4, 1)
    #R of the species. Eventually can vary through time according to env condition
  b <- .8#runif(1, .1, .9)
   #Density dependence on the species on itself. Higher = less DD
  
  #Create the model
  for (t in 1:(time-1)) {
    #calculate population sizes for species 1
    N1[t+1] <- r + b*N1[t]
    
    #if(is.nan(N1[t+1])) {
      #N1[t+1] <- 0
    #}
    
    #if(N1[t+1] < 1) {
      #N1[t+1] <- 0
    #}
    plot(N1, type = "p")
  }
#Always goes to extinction within 10 timesteps


  
  
  
  
########
#Multiple species
########
#runs <- 20 #number of runs
time <- 30
species <- 5

#Set parameters
r <- runif(species, 0, 1)
  #Vector of r for each species. Eventually can vary through time according to env condition?
#b <- runif(species, .1, 1)

#Set b to a specific number
b <- rep(0.6, species)   #Use consistent 
b <- diag(b)
b <- ifelse(0, runif(species*species, -.25, .25), b)
  ##Density dependence on the species on itself. Higher = less DD

#Set state variables
N <- matrix(NA, nrow = time, ncol = species) 
#Starting abundance per species
#N[1,] <- log(5)
N[1,] <- -3

#set.seed(555)
#Create the model
for (t in 1:(time-1)) {
  
  #E <- rnorm(species, mean = 0, sd = 1)
  
  #for(s in 1:species){
    #calculate population sizes for all species
    N[t+1,] <- r + b%*%N[t,] #+ E
  
    #if(is.nan(N[t+1,s])) {
      #N[t+1,s] <- 0
    #}
 # }
  
  Nplot <- pivot_longer(as.data.frame(N), cols = 1:species, names_to = "Species", 
                        values_to = "Abundance") %>% 
    mutate(Time = rep(1:time, each = species))

}

ggplot(Nplot, aes(x = Time, y = Abundance, color = Species)) +
  geom_point() +
  geom_line()
