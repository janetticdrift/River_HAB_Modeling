#Historical Predictions, 2022
library(ggplot2)

#Pull out community abundances and demographics 
x <- rstan::extract(fit.m4)
abundances <- x[["n"]][c(1:10),,1]
alphas <- x[["Alpha"]][c(1:10),]
betas <- as.array(x[["Beta"]])[c(1:10),,] 
sigmas <- x[["sigma_p"]][c(1:10),]

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(time, 4, runs))

#Pull out environmental effects
Ntheta <- x[["Ntheta"]][c(1:10),]
nitrate <- stand_nut$nitrate_mg_N_L[1:time]

Ptheta <- x[["Ptheta"]][c(1:10),]
phos <- stand_nut$oPhos_ug_P_L[1:time]

Atheta <- x[["Atheta"]][c(1:10),]
amon <- stand_nut$ammonium_mg_N_L[1:time]

Dtheta <- x[["Dtheta"]][c(1:10),]
dis <- discharge$stand_discharge[1:time]

Ttheta <- x[["Ttheta"]][c(1:10),]
temp <- stand_nut$temp_C[1:time]


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
  
  
  for(t in 2:time){
      n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z] + nTheta*nitrate[t-1]+
                                 pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] + 
                                 tTheta*temp[t-1], Sigma = sigma)
  }
}

sims <- as.data.frame(apply(n, c(1,2), mean)) %>% 
  mutate(across(1:4, exp)) %>%
  dplyr::rename(green_algae = V1, microcoleus = V2,
                anabaena_cylindrospermum = V3,
                other_nfixers = V4) %>% 
  mutate(time = 1:13) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance") %>% 
  mutate(real_week = time + 25, year = 2022) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



ggplot(sims, aes(x = model_date, y = Abundance, colour = Species)) +
  geom_line() +
  scale_y_continuous(breaks=c(seq(0,100,5)))
  
