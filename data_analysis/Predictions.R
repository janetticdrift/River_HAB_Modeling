#Historical Predictions, 2022

#Pull out starting abundances 
x <- rstan::extract(fit.m4)
abundances <- x[["n"]][c(1:100),,1]
alphas <- x[["Alpha"]][c(1:100),]
betas <- as.array(x[["Beta"]]) %>% 
  head(100)

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022

for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[z,]
  Beta <- betas[z,,]
  n <- abundances[z,]
  
  for(t in 2:time){ 
      n[,t] ~ multi_normal(Alpha + Beta*n[,t-1], ID)
  }
}
