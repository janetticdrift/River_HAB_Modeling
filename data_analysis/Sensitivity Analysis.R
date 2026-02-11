###########################
#Sensitivity Analyses: Finding Equilibrium Abundance
###########################

#Read in model output
  #Model used includes all biotic and abiotic variables
x <- rstan::extract(fit.m4) #m4 = all vars

#Pull out species interaction effects
betas <- as.array(x[["Beta"]])[,,]

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
env_range <- seq(-3, 3, by = 0.5)
env_names <- names(theta_list)

w_star <- array(        # iteration, env_var, perturbation, species
  NA, 
  dim = c(15000, length(env_names), length(env_range), 4),
  dimnames = list(
    iteration = NULL,
    env = env_names,
    env_peturb = env_range,
    species = paste0("sp", 1:4)
  )
)

#Run model
for(z in 1:15000){
  
  Beta <- betas[z,,]
  
  #Include species interactions
  M_inverse <- solve(ID - Beta)    #solve() takes the inverse of matrices
  
  #Exclude species interactions
  # B_intra <- Beta
  # B_intra[row(B_intra) != col(B_intra)] <- 0    #any time row doesn't equal column, set to 0
  # M_inv_intra <- solve(ID - B_intra)
  
  for(v in seq_along(env_names)){     #Use seq_along for iterating along characters
    
    theta_matrix <- theta_list[[v]]   #Pull out current env var of interest
    theta_vec <- theta_matrix[z,]     #Pull out theta vector of z iteration
    
    for(e in seq_along(env_range)){
      
      A <- theta_vec * env_range[e]
      w_star[z, v, e, ] <- as.vector(M_inverse %*% A) #Switch between M_inverse and M_inv_intra
      
    }
  }
}

#Clean output for plotting
eq_abund <- as.data.frame.table(w_star, responseName = "w_star") %>%
  mutate(
    env_peturb = as.numeric(as.character(env_peturb)),
    iteration = as.integer(iteration),
    species = factor(species), # Green algae, Microcoleus, Anabaena, Other N-fixers
    env = factor(env)
  )

#Plot outputs
  #Species2: Microcoleus
ggplot(subset(eq_abund, species %in% "sp2"), aes(x = factor(env_peturb), y = w_star)) +
  facet_wrap(~env, scales = "free") +
  geom_boxplot(outliers = T)


