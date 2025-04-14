#Historical Predictions, 2022
library(ggplot2)

#Pull out starting abundances 
x <- rstan::extract(fit.m4)
abundances <- x[["n"]][c(1:10),,1]
alphas <- x[["Alpha"]][c(1:10),]
betas <- as.array(x[["Beta"]])[c(1:10),,] #%>% 
  #head(10)

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(time, 4, runs))

for(z in 1:runs){
  #Set parameters
  Alpha <- alphas[1,]
  Beta <- betas[1,,]
  n[1,,z] <- abundances[1,]
  sigma <- diag(rnorm(4, mean = 0, sd = 1))
  
  for(t in 2:time){
      n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta*n[t-1,,z], Sigma = sigma[z])
  }
}

n[2,,] <- MASS::mvrnorm(n = 15, mu = matrix(c(1:4, 5:8), nrow = 4, ncol = 4), Sigma = sigma)
n[2,,1] <- Rfast::mvnorm.mle(Alpha + Beta*n[1,,1])

print(n)

sims <- as.data.frame(apply(n, c(1,2), mean)) %>% 
  mutate(across(1:4, exp)) %>%
  dplyr::rename(green_algae = V1, microcoleus = V2,
                anabaena_cylindrospermum = V3,
                other_nfixers = V4) %>% 
  mutate(time = 1:13) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")



ggplot(sims, aes(x = time, y = Abundance, colour = Species)) +
  geom_line()
  
