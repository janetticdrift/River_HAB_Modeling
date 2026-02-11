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

#Create identity matrix
ID <- diag(1, nrow = 4, ncol = 4) #4 is the number of species in the river-wide model

#Find equilibrium abundances
#Equi_abund_W = A / (Identity matrix - B)

for(z in 1:runs){
  #Set parameters
  Beta <- betas[z,,]
  
  #Pull env covariates
  nTheta <- Ntheta[z,]
  pTheta <- Ptheta[z,]
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  #Model
  w[,,z] <- nTheta / (ID - Beta)
    n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z],
                             Sigma = sigma)
}