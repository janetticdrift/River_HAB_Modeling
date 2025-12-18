#Historical Predictions, 2022
library(ggplot2)
library(ggpubr)
library(tidyverse)

#Pull out community abundances and demographics 
x <- rstan::extract(fit.m1) #m1 = averaged reaches for Microcoleus
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- as.array(x[["Beta"]])[,,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #13 weeks in 2022, 13 in 2023, 15 in 2024
n <- array(NA, dim = c(time, 9, runs)) #9 is number of species

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

Rtheta <- x[["Rtheta"]][,]
rad <- swradiation$stand_rad[1:time]

#Run predictions
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
    #Include env drivers
    n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    # #Remove env drivers
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z],
    #                          Sigma = sigma)
  }
}

matsims2022 <- as.data.frame(apply(n, c(1,2), mean)) %>% 
  dplyr::mutate(across(1:9, exp)) %>%
  dplyr::rename(anabaena_and_cylindrospermum = V1, 
                e_diatoms = V2,
                geitlerinema = V3,
                green_algae = V4, 
                microcoleus = V5,
                non_e_diatoms = V6,
                nostoc = V7,
                other_coccoids = V8,
                rare = V9) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:9, names_to = "Species", values_to = "Abundance") %>% 
  mutate(real_week = time + 25, year = 2022) %>% #HEY THIS CHANGED FROM 24 TO 25
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))

#Plot
mat_p22 <- ggplot(matsims2022, aes(x = model_date, y = Abundance, colour = Species)) +
  geom_line(size = 1.5) +
  geom_line(data=mat_params2[mat_params2$year %in% "2022", ], 
            aes(x = model_date, y = mean, colour = Species),
            linewidth = 4, alpha = .25) +
  # geom_ribbon(aes(ymin = `CIlower`, ymax = `CIupper`,
  #                 fill = Species), alpha = 0.3) +
  scale_y_continuous(breaks=c(seq(0,150,5))) +
  #coord_cartesian(ylim = c(0,70)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "2022 Predictions") +
  scale_color_manual(labels = c("Anabaena", "Epithemia", "Geitlerinema", 
                                "Green Algae", "Microcoleus","Non-Epithemia",
                                "Nostoc","Other Coccoids", "Rare"), values = c("darkcyan", "darkgreen", 
                                                               "darkorange", "darkolivegreen4", "brown", 
                                                               "red", "darkviolet","navy","goldenrod"))

#####################################
#Historical Predictions, 2023
#####################################
#Pull out community abundances and demographics 
abundances <- x[["n"]][,,14] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2023
n <- array(NA, dim = c(time, 9, runs))

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[14:(13+time)]

phos <- stand_nut$oPhos_ug_P_L[14:(13+time)]

amon <- stand_nut$ammonium_mg_N_L[14:(13+time)]

dis <- discharge$stand_discharge[14:(13+time)]

temp <- stand_nut$temp_C[14:(13+time)]

cond <- stand_nut$cond_uS_cm[14:(13+time)]

rad <- swradiation$stand_rad[14:(13+time)]


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
    #Everything included
    n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    # #Remove env drivers
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z],
    #                          Sigma = sigma)
    
    
  }
}

matsims2023 <- as.data.frame(apply(n, c(1,2), mean)) %>% 
  dplyr::mutate(across(1:9, exp)) %>%
  dplyr::rename(anabaena_and_cylindrospermum = V1, 
                e_diatoms = V2,
                geitlerinema = V3,
                green_algae = V4, 
                microcoleus = V5,
                non_e_diatoms = V6,
                nostoc = V7,
                other_coccoids = V8,
                rare = V9) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:9, names_to = "Species", values_to = "Abundance") %>% 
  mutate(real_week = time + 26, year = 2023) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



mat_p23 <- ggplot(matsims2023, aes(x = model_date, y = Abundance, colour = Species)) +
  geom_line(size = 1.5) +
  geom_line(data=mat_params2[mat_params2$year %in% "2023", ], 
            aes(x = model_date, y = mean, colour = Species),
            linewidth = 4, alpha = .25) +
  scale_y_continuous(breaks=c(seq(0,100,5))) +
  #coord_cartesian(ylim = c(0,50)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "2023 Predictions") +
  scale_color_manual(labels = c("Anabaena", "Epithemia", "Geitlerinema", 
                                "Green Algae", "Microcoleus","Non-Epithemia",
                                "Nostoc","Other Coccoids", "Rare"), values = c("darkcyan", "darkgreen", 
                                                                               "darkorange", "darkolivegreen4", "brown", 
                                                                               "red", "darkviolet","navy","goldenrod"))

#####################################
#Historical Predictions, 2024
#####################################
#Pull out community abundances and demographics 
abundances <- x[["n"]][,,27] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(time, 9, runs))

#Pull out environmental effects
nitrate <- stand_nut$nitrate_mg_N_L[27:(26+time)]

phos <- stand_nut$oPhos_ug_P_L[27:(26+time)]

amon <- stand_nut$ammonium_mg_N_L[27:(26+time)]

dis <- discharge$stand_discharge[27:(26+time)]

temp <- stand_nut$temp_C[27:(26+time)]

cond <- stand_nut$cond_uS_cm[27:(26+time)]

rad <- swradiation$stand_rad[27:(26+time)]


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
    
    #Everything included
    n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    # #Remove env drivers
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z],
    #                          Sigma = sigma)
    
  }
}

matsims2024 <- as.data.frame(apply(n, c(1,2), mean)) %>% 
  dplyr::mutate(across(1:9, exp)) %>%
  dplyr::rename(anabaena_and_cylindrospermum = V1, 
                e_diatoms = V2,
                geitlerinema = V3,
                green_algae = V4, 
                microcoleus = V5,
                non_e_diatoms = V6,
                nostoc = V7,
                other_coccoids = V8,
                rare = V9) %>%  
  mutate(time = 1:time) %>% 
  pivot_longer(cols = 1:9, names_to = "Species", values_to = "Abundance") %>% 
  mutate(real_week = time + 26, year = 2024) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))



mat_p24 <- ggplot(matsims2024, aes(x = model_date, y = Abundance, colour = Species)) +
  #geom_point(size = 3)+
  geom_line(size = 1.5) +
  geom_line(data=mat_params2[mat_params2$year %in% "2024", ], 
            aes(x = model_date, y = mean, colour = Species),
            linewidth = 4, alpha = .25) +
  scale_y_continuous(breaks=c(seq(0,100,5))) +
  #coord_cartesian(ylim = c(0,45)) +
  labs(x = "Date", y = "Percent Cover (%)", title = "2024 Predictions") +
  scale_color_manual(labels = c("Anabaena", "Epithemia", "Geitlerinema", 
                                "Green Algae", "Microcoleus","Non-Epithemia",
                                "Nostoc","Other Coccoids", "Rare"), values = c("darkcyan", "darkgreen", 
                                                                               "darkorange", "darkolivegreen4", "brown", 
                                                                               "red", "darkviolet","navy","goldenrod"))



ggarrange(
  mat_p22, mat_p23, mat_p24, labels = c("A", "B", "C"), ncol = 3,
  common.legend = TRUE, legend = "bottom"
)




