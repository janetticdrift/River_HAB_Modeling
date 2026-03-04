###########################
#Sensitivity Analyses: Finding Equilibrium Abundance
###########################

#Read in model output
  #Model used includes all biotic and abiotic variables
x <- rstan::extract(fit.m4) #m4 = all vars

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
env_range <- seq(-3, 3, length.out = 11)   #Range of env effect, in standard deviation units
env_names <- names(theta_list)      #Names of env variables

w_star_include <- array(        # iteration, env_var, perturbation, species
  NA, 
  dim = c(15000, length(env_names), length(env_range), 4),
  dimnames = list(
    iteration = NULL,
    env = env_names,
    env_peturb = env_range,
    species = 1:4
  )
)

#Run model
for(z in 1:15000){
  
  Alpha <- alphas[z,]   #Current iteration of Alpha: baseline growth
  Beta <- betas[z,,]    #Current iteration of Beta: species interactions
  
  #Include species interactions
  M_inverse <- solve(ID - Beta)    #solve() takes the inverse of matrices
  
  for(v in seq_along(env_names)){     #Use seq_along for iterating along characters
    
    theta_matrix <- theta_list[[v]]   #Pull out current env var of interest
    theta_vec <- theta_matrix[z,]     #Pull out theta vector of z iteration
    
    for(e in seq_along(env_range)){
      
      A <- (theta_vec * env_range[e]) + Alpha     #Multiply coef by std dev, then add intercept
      w_star_include[z, v, e, ] <- as.vector(M_inverse %*% A) #Switch between M_inverse and M_inv_intra
      
    }
  }
}

#Clean output for plotting
eq_abund_include <- as.data.frame.table(w_star_include, responseName = "w_star") %>%
  mutate(
    env_peturb = as.numeric(as.character(env_peturb)),
    iteration = as.integer(iteration),
    species = factor(species), # Green algae, Microcoleus, Anabaena, Other N-fixers
    env = factor(env)
    ) %>% 
  select(!iteration) %>% 
  mutate(Interactions = "include") %>% 
  mutate(w_star = exp(w_star)) %>%  # Back transform w_star from logged values
                               #About 200 values were removed for being too large: exp 
                               #returns Inf on values over 710
  filter(is.finite(w_star)) #Remove those Inf values

#Read in model output-------------------------------------------------------------
#Model used includes only abiotic variables
x <- rstan::extract(fit.m6) #m6 = all vars

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
  dim = c(15000, length(env_names), length(env_range), 4),
  dimnames = list(
    iteration = NULL,
    env = env_names,
    env_peturb = env_range,
    species = 1:4
  )
)

#Run model
for(z in 1:15000){
  
  Alpha <- alphas[z,]   #Current iteration of Alpha: baseline growth
  Beta <- betas[z,]    #Current iteration of Beta: species interactions

  #Exclude species interactions
  M_inverse <- solve(ID - B_intra)
  
  for(v in seq_along(env_names)){     #Use seq_along for iterating along characters
    
    theta_matrix <- theta_list[[v]]   #Pull out current env var of interest
    theta_vec <- theta_matrix[z,]     #Pull out theta vector of zth iteration
    
    for(e in seq_along(env_range)){
      
      A <- (theta_vec * env_range[e]) + Alpha     #Multiply coef by env_peturb, then add intercept
      w_star_exclude[z, v, e, ] <- as.vector(M_inverse %*% A) #Switch between M_inverse and M_inv_intra
      
    }
  }
}

#Clean output for plotting
eq_abund_exclude <- as.data.frame.table(w_star_exclude, responseName = "w_star") %>%
  mutate(
    env_peturb = as.numeric(as.character(env_peturb)),
    iteration = as.integer(iteration),
    species = factor(species), # Green algae, Microcoleus, Anabaena, Other N-fixers
    env = factor(env)) %>% 
  select(!iteration) %>% 
  mutate(Interactions = "exclude") %>% 
  mutate(w_star = exp(w_star)) #Model excluding sp interaction did not return values too large

#Bind together predictions where species were included or excluded in predictions
  # Function for ID'ing outliers
outlier <- function(x) {
  # Calculate first and third quartiles
  Q1 <- quantile(x, probs = .25, na.rm = TRUE)
  Q3 <- quantile(x, probs = .75, na.rm = TRUE)
  # Calculate the Interquartile Range (IQR)
  IQR <- Q3 - Q1
  # Define lower and upper bounds
  lower_bound <- Q1 - (1.5 * IQR)
  upper_bound <- Q3 + (1.5 * IQR)
  # Return vector indicating which values are NOT outliers
  x >= lower_bound & x <= upper_bound
}
eq_abund <- rbind(eq_abund_include, eq_abund_exclude) %>% 
  group_by(Interactions) %>% 
  filter(outlier(w_star))



#Plot outputs---------------------------------------------------------------
  #Boxplots
env_labels <- c(
  'N' = 'Nitrate',
  'P' = 'Phosphate',
  'A' = 'Ammonium',
  'D' = 'Discharge',
  'Tt' = 'Temperature',
  'C' = 'Conductivity',
  'R' = 'Light'
)
  #Species2: Microcoleus
ggplot(subset(eq_abund, species %in% "2" & Interactions %in% "include"), aes(x = factor(env_peturb), y = w_star, fill = env_peturb)) +
  facet_wrap(~env, labeller = labeller(env = env_labels)) +
  geom_boxplot(outliers = F) + 
  labs(x = "Standard Deviations", y = "Equilibrium Abundance") +
  scale_x_discrete(breaks = c(-3, 0, 3)) +
  labs(title = "Microcoleus") +
  theme_classic() +
  scale_fill_viridis_c(option="magma", begin = 0.45, end = .80)+
  theme(legend.position = "none")

#Species2: Anabaena
ggplot(subset(eq_abund, species %in% "3" & Interactions %in% "include"), aes(x = factor(env_peturb), y = w_star, fill = env_peturb)) +
  facet_wrap(~env, labeller = labeller(env = env_labels)) +
  geom_boxplot(outliers = F) + 
  labs(x = "Standard Deviations", y = "Equilibrium Abundance") +
  scale_x_discrete(breaks = c(-3, 0, 3)) +
  labs(title = "Anabaena") +
  theme_classic() +
  scale_fill_viridis_c(option="viridis", begin = 0.5, end = .90)+
  theme(legend.position = "none")

  #Regression plot 
  #Species2: Microcoleus
ggplot(subset(eq_abund %>% slice_sample(n = 10000), species %in% "2"), aes(x = env_peturb, y = w_star, color = Interactions)) +
  facet_wrap(~env, labeller = labeller(env = env_labels)) +
  geom_smooth(method = "loess", se = F) +
  labs(x = "Standard Deviations", y = "Equilibrium Abundance", title = "Microcoleus") +
  guides(color = guide_legend(reverse = TRUE)) +
  scale_x_continuous(breaks = c(-3, 0, 3)) +
  theme_classic() +
  scale_color_manual(values = c("#E69F00", "#56B4E9"))

#Species2: Anabaena
ggplot(subset(eq_abund %>% slice_sample(n = 10000), species %in% "3"), aes(x = env_peturb, y = w_star, color = Interactions)) +
  facet_wrap(~env, labeller = labeller(env = env_labels)) +
  geom_smooth(method = "loess") +
  labs(x = "Standard Deviations", y = "Equilibrium Abundance", title = "Anabaena") +
  guides(color = guide_legend(reverse = TRUE)) +
  scale_x_continuous(breaks = c(-3, 0, 3)) +
  theme_classic() +
  scale_color_manual(values = c("#E69F00", "#56B4E9"))


any(is.infinite(eq_abund_exclude$w_star))
any(is.na(eq_abund$w_star))

