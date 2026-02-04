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

#Find equilibrium abundances
#Equi_abund_W = A / (Identity matrix - B)

for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,,]
  n[1,,z] <- abundances[z,]
  sigma <- diag(sigmas[z,])
  
  #Pull env covariates
  nTheta <- Ntheta[z,]
  pTheta <- Ptheta[z,]
  aTheta <- Atheta[z,]
  dTheta <- Dtheta[z,]
  tTheta <- Ttheta[z,]
  cTheta <- Ctheta[z,]
  rTheta <- Rtheta[z,]
  
  
  for(t in 2:time){
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z], Sigma = sigma)
    
    #Everything included
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z] + nTheta*nitrate[t-1]+
    #                            pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
    #                            tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
    #                          Sigma = sigma)
    # #Remove env drivers
    n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z],
                             Sigma = sigma)
  }
}