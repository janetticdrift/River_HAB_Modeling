#Historical Predictions, 2022
library(ggplot2)

#Pull out community abundances and demographics 
x <- rstan::extract(fit.m4)
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- as.array(x[["Beta"]])[,,] 
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2022
n <- array(NA, dim = c(time, 4, runs))

#Pull out environmental effects
Ntheta <- x[["Ntheta"]][,]
nitrate <- stand_nut$nitrate_mg_N_L[1:time]

Ptheta <- x[["Ptheta"]][,]
phos <- stand_nut$oPhos_ug_P_L[1:time]

Atheta <- x[["Atheta"]][,]
amon <- stand_nut$ammonium_mg_N_L[1:time]

Dtheta <- x[["Dtheta"]][,]
dis <- discharge$stand_discharge[1:time]

Ttheta <- x[["Ttheta"]][,]
temp <- stand_nut$temp_C[1:time]

Ctheta <- x[["Ctheta"]][,]
cond <- stand_nut$cond_uS_cm[1:time]


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
  
  
  for(t in 2:time){
      # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z], Sigma = sigma)
      
      n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z] + nTheta*nitrate[t-1]+
                                 pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                                 tTheta*temp[t-1] + cTheta*cond[t-1], Sigma = sigma)
  }
}

#Create dataframe
sims2022mean <- as.data.frame(apply(n, c(1,2), mean)) %>% 
  mutate(across(1:4, exp)) %>%
  dplyr::rename(green_algae = V1, microcoleus = V2,
                anabaena_cylindrospermum = V3,
                other_nfixers = V4) %>% 
  mutate(time = 1:time) %>%
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance")
sims2022lquant <- as.data.frame(apply(n, c(1,2), quantile, probs = 0.025)) %>% 
  mutate(across(1:4, exp)) %>%
  dplyr::rename(green_algae = V1, microcoleus = V2,
                anabaena_cylindrospermum = V3,
                other_nfixers = V4) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIlower")
sims2022uquant <- as.data.frame(apply(n, c(1,2), quantile, probs = 0.975)) %>% 
  mutate(across(1:4, exp)) %>%
  dplyr::rename(green_algae = V1, microcoleus = V2,
                anabaena_cylindrospermum = V3,
                other_nfixers = V4) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "CIupper")

sims2022 <- left_join(sims2022mean, sims2022lquant, by=c("Species", "time")) %>%
                        left_join(., sims2022uquant, by=c("Species", "time")) %>% 
  mutate(real_week = time + 25, year = 2022) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))


#Plot
p22 <- ggplot(sims2022, aes(x = model_date, y = Abundance, colour = Species)) +
  geom_line(size = 1) +
  # geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`,
  #                 fill = Species), alpha = 0.3) +
  scale_y_continuous(breaks=c(seq(0,150,5))) +
  coord_cartesian(ylim = c(0,70)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "2022 Predictions") 
  



#####################################
#Historical Predictions, 2023
#####################################
#Pull out community abundances and demographics 
abundances <- x[["n"]][,,14] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(time, 4, runs))

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[1:time]

phos <- stand_nut$oPhos_ug_P_L[1:time]

amon <- stand_nut$ammonium_mg_N_L[1:time]

dis <- discharge$stand_discharge[1:time]

temp <- stand_nut$temp_C[1:time]

cond <- stand_nut$cond_uS_cm[1:time]


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
  
  
  for(t in 2:time){
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z], Sigma = sigma)
    
    n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1], Sigma = sigma)
  }
}

sims2023 <- as.data.frame(apply(n, c(1,2), mean)) %>% 
  mutate(across(1:4, exp)) %>%
  dplyr::rename(green_algae = V1, microcoleus = V2,
                anabaena_cylindrospermum = V3,
                other_nfixers = V4) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance") %>% 
  mutate(real_week = time + 24, year = 2023) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



p23 <- ggplot(sims2023, aes(x = model_date, y = Abundance, colour = Species)) +
  geom_line(size = 1) +
  scale_y_continuous(breaks=c(seq(0,100,5))) +
  coord_cartesian(ylim = c(0,40)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "2023 Predictions") 



#####################################
#Historical Predictions, 2024
#####################################
#Pull out community abundances and demographics 
abundances <- x[["n"]][,,29] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 17 #number of weeks in 2023
n <- array(NA, dim = c(time, 4, runs))

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[1:time]

phos <- stand_nut$oPhos_ug_P_L[1:time]

amon <- stand_nut$ammonium_mg_N_L[1:time]

dis <- discharge$stand_discharge[1:time]

temp <- stand_nut$temp_C[1:time]

cond <- stand_nut$cond_uS_cm[1:time]


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
  
  
  for(t in 2:time){
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z], Sigma = sigma)
    
    n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1], Sigma = sigma)
  }
}

sims2024 <- as.data.frame(apply(n, c(1,2), mean)) %>% 
  mutate(across(1:4, exp)) %>%
  dplyr::rename(green_algae = V1, microcoleus = V2,
                anabaena_cylindrospermum = V3,
                other_nfixers = V4) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:4, names_to = "Species", values_to = "Abundance") %>% 
  mutate(real_week = time + 24, year = 2024) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



p24 <- ggplot(sims2024, aes(x = model_date, y = Abundance, colour = Species)) +
  geom_line(size = 1) +
  scale_y_continuous(breaks=c(seq(0,100,5))) +
  coord_cartesian(ylim = c(0,35)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "2024 Predictions") 


ggarrange(
  p22, p23, p24, labels = c("A", "B", "C"), ncol = 3,
  common.legend = TRUE, legend = "bottom"
)

