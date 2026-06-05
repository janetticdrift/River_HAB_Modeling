###########################
#Sensitivity Analyses: Finding Equilibrium Abundance
###########################
library(tidyverse)
library(here)
library(ggplot2)
library(patchwork)

#Read in model output
allvar <- readRDS(here::here("data/Riverwide_AllVar_predictions.rds")) #Model used includes all biotic and abiotic variables
abiotic <- readRDS(here::here("data/Riverwide_Abiotic_predictions.rds")) #Model used includes abiotic variables


###############
#Calculate w* for model including biotic interactions
x <- allvar

#Pull out species demographics
betas <- as.array(x[["Beta"]])[,,]
alphas <- x[["Alpha"]][,]

#Pull out environmental effects
Ntheta <- x[["Ntheta"]][,]
Ptheta <- x[["Ptheta"]][,]
Atheta <- x[["Atheta"]][,]
Dtheta <- x[["Dtheta"]][,]
Ttheta <- x[["Ttheta"]][,]
Ctheta <- x[["Ctheta"]][,]
Rtheta <- x[["Rtheta"]][,]
#Combine env theta effects into a list
theta_list <- list(
  N = Ntheta,
  P = Ptheta,
  A = Atheta,
  D = Dtheta,
  Tt = Ttheta,
  C = Ctheta,
  R = Rtheta
)

#Create identity matrix
ID <- diag(1, nrow = 4, ncol = 4) # 4 is the number of species in the river-wide model

#Set conditions and output storage
env_range <- seq(-2, 2, length.out = 11)   #Range of env effect, in standard deviation units
env_names <- names(theta_list)      #Names of env variables

w_star_include <- array(        # iteration, env_var, perturbation, species
  NA, 
  dim = c(9000, length(env_names), length(env_range), 4),
  dimnames = list(
    iteration = NULL,
    env = env_names,
    env_peturb = env_range,
    species = 1:4
  )
)

#Run model
for(z in 1:9000){
  
  Alpha <- alphas[z,]   #Current iteration of Alpha: baseline growth
  Beta <- betas[z,,]    #Current iteration of Beta: species interactions
  
  #Include species interactions
  M_inverse <- solve(ID - Beta)    #solve() takes the inverse of matrices
  
  for(v in seq_along(env_names)){     #Use seq_along for iterating along character vectors
    
    theta_matrix <- theta_list[[v]]   #Pull out current env var of interest
    theta_vec <- theta_matrix[z,]     #Pull out env theta coefficients of z iteration
    
    for(e in seq_along(env_range)){
      
      A <- (theta_vec * env_range[e]) + Alpha  #Multiply env coef by perturbation amount, then add intercept
      w_star_include[z, v, e, ] <- as.vector(M_inverse %*% A) #Calculate w*
      
    }
  }
}

#Clean output for plotting
eq_abund_include <- as.data.frame.table(w_star_include, responseName = "w_star") %>%
  dplyr::mutate(
    env_peturb = as.numeric(as.character(env_peturb)),
    iteration = as.integer(iteration),
    species = factor(species), # Green algae, Microcoleus, Anabaena, Other N-fixers
    env = factor(env)
    ) %>% 
  dplyr::select(!iteration) %>% 
  dplyr::mutate(Interactions = "include") %>% #Add a category column for including/excluding species interactions
  dplyr::mutate(w_star = exp(w_star)) %>%  # Back transform w_star from logged values
                               # Values over 710 result in -Inf values when exp'ed for being too large
  dplyr::filter(is.finite(w_star)) #Remove those Inf values

#Read in model output-------------------------------------------------------------
#Model used includes only abiotic variables

x <- abiotic 

#Pull out species demographics
betas <- as.array(x[["Beta"]])[,]
alphas <- x[["Alpha"]][,]

#Pull out environmental effects
Ntheta <- x[["Ntheta"]][,]
Ptheta <- x[["Ptheta"]][,]
Atheta <- x[["Atheta"]][,]
Dtheta <- x[["Dtheta"]][,]
Ttheta <- x[["Ttheta"]][,]
Ctheta <- x[["Ctheta"]][,]
Rtheta <- x[["Rtheta"]][,]
#Combine env theta effects into a list
theta_list <- list(
  N = Ntheta,
  P = Ptheta,
  A = Atheta,
  D = Dtheta,
  Tt = Ttheta,
  C = Ctheta,
  R = Rtheta
)

w_star_exclude <- array(        # iteration, env_var, perturbation, species
  NA, 
  dim = c(9000, length(env_names), length(env_range), 4),
  dimnames = list(
    iteration = NULL,
    env = env_names,
    env_peturb = env_range,
    species = 1:4
  )
)

#Run model
for(z in 1:9000){
  
  Alpha <- alphas[z,]   #Current iteration of Alpha: baseline growth
  Beta_vec <- betas[z,]    #Current iteration of Beta: species interactions
  Beta <- diag(Beta_vec)  #Convert interactions into matrix form
  
  #Exclude species interactions
  M_inverse <- solve(ID - Beta)
  
  for(v in seq_along(env_names)){     #Use seq_along for iterating along characters
    
    theta_matrix <- theta_list[[v]]   #Pull out current env var of interest
    theta_vec <- theta_matrix[z,]     #Pull out theta vector of zth iteration
    
    for(e in seq_along(env_range)){
      
      A <- (theta_vec * env_range[e]) + Alpha     #Multiply coef by env_peturb, then add intercept
      w_star_exclude[z, v, e, ] <- as.vector(M_inverse %*% A) 
      
    }
  }
}

#Clean output for plotting
eq_abund_exclude <- as.data.frame.table(w_star_exclude, responseName = "w_star") %>%
  dplyr::mutate(
    env_peturb = as.numeric(as.character(env_peturb)),
    iteration = as.integer(iteration),
    species = factor(species), 
    env = factor(env)) %>% 
  dplyr::select(!iteration) %>% 
  dplyr::mutate(Interactions = "exclude") %>% 
  dplyr::mutate(w_star = exp(w_star)) %>% 
  dplyr::filter(is.finite(w_star))


#Bind together predictions where species were included or excluded in predictions
eq_abund <- rbind(eq_abund_include, eq_abund_exclude) %>% 
  group_by(Interactions)


#Checking for errors
any(is.infinite(eq_abund_exclude$w_star))
any(is.na(eq_abund$w_star))



#Plot outputs---------------------------------------------------------------

#Create an object that renames the environmental variables from their abbreviations
env_labels <- c(
  'N' = 'Nitrate',
  'P' = 'Phosphate',
  'A' = 'Ammonium',
  'D' = 'Discharge',
  'Tt' = 'Temperature',
  'C' = 'Conductivity',
  'R' = 'Light'
)

#########  
#Boxplots
######### 
  #Species2: Microcoleus
ggplot(subset(eq_abund, species %in% "2" & Interactions %in% "include"), 
       aes(x = factor(env_peturb), y = w_star, fill = env_peturb)) +
  facet_wrap(~env, scales = "free_y", labeller = labeller(env = env_labels)) +
  geom_boxplot(outliers = F) + 
  labs(x = "Standard Deviations", y = "Equilibrium Abundance") +
  scale_x_discrete(breaks = c(-2, 0, 2)) +
  labs(title = "Microcoleus") +
  theme_classic() +
  scale_fill_viridis_c(option="magma", begin = 0.45, end = .80)+
  theme(legend.position = "none")

#Species2: Anabaena
ggplot(subset(eq_abund, species %in% "3" & Interactions %in% "include"), aes(x = factor(env_peturb), y = w_star, fill = env_peturb)) +
  facet_wrap(~env, scales = "free_y", labeller = labeller(env = env_labels)) +
  geom_boxplot(outliers = F) + 
  labs(x = "Standard Deviations", y = "Equilibrium Abundance") +
  scale_x_discrete(breaks = c(-2, 0, 2)) +
  labs(title = "Anabaena") +
  theme_classic() +
  scale_fill_viridis_c(option="viridis", begin = 0.5, end = .90)+
  theme(legend.position = "none")

######### 
#Median plot 
######### 

eq_abund_medians <- eq_abund %>%
  group_by(env, env_peturb, species, Interactions) %>%
  dplyr::summarise(med_star = median(w_star, na.rm = TRUE))

  #Species2: Microcoleus
ggplot(subset(eq_abund_medians, species %in% "2"), aes(x = env_peturb, y = med_star, color = Interactions)) +
  facet_wrap(~env, labeller = labeller(env = env_labels)) +
  geom_line() +
  guides(color = guide_legend(reverse = TRUE)) +
  labs(x = "Standard Deviations", y = "Equilibrium Abundance", title = "Microcoleus") +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  theme_classic()

#Species2: Anabaena
ggplot(subset(eq_abund_medians, species %in% "3"), aes(x = env_peturb, y = med_star, color = Interactions)) +
  facet_wrap(~env, labeller = labeller(env = env_labels)) +
  geom_line() +
  guides(color = guide_legend(reverse = TRUE)) +
  labs(x = "Standard Deviations", y = "Equilibrium Abundance", title = "Anabaena") +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  theme_classic()

