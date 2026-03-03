#Historical Predictions, 2022
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(abind)

#files mat_params2 and mat_params2_groups are in WithinMatModel_vs_real
  

#Pull out community abundances and demographics 
x <- rstan::extract(fit.m1) #m1 = averaged reaches for Microcoleus
abundances <- x[["n"]][,,1] #iterations, species #, time
alphas <- x[["Alpha"]][,]
betas <- as.array(x[["Beta"]])[,,]
sigmas <- x[["sigma_p"]][,]

#inputs
runs <- nrow(abundances)
time <- 13 #13 weeks in 2022, 13 in 2023, 15 in 2024
n <- array(NA, dim = c(runs, 11, time)) #11 is number of species

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
  n[z,,1] <- abundances[z,]
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
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    # #Remove env drivers
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z],
    #                          Sigma = sigma)
  }
}

pred2022 <- n

matsims2022 <- as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:11, exp)) %>%
  dplyr::rename(Anabaena = V1, 
                'Epithemia Diatoms' = V2,
                Geitlerinema = V3,
                'Green Algae' = V4, 
                Leptolyngbya = V5,
                Microcoleus = V6,
                'Non-Epithemia Diatoms' = V7,
                Nostoc = V8,
                Oscillatoria = V9,
                'Other Coccoids' = V10,
                Rare = V11) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = c(1:11), names_to = "Species", values_to = "Abundance") %>% #cols are 1:9 if species are ungrouped, 1:3, 5:6 if grouped
  mutate(real_week = time + 25, year = 2022) %>% #HEY THIS CHANGED FROM 24 TO 25
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))

#Create a color palette
mycols <- c("brown", "chartreuse3", "goldenrod", "darkolivegreen4", 
            "violetred", "darkcyan", "mediumpurple3","khaki1", "antiquewhite3",
            "darkorange", "lightblue1")
mypal <- palette(mycols)
names(mypal) = c("Anabaena", "Epithemia Diatoms", "Geitlerinema", 
                 "Green Algae", "Leptolyngbya", "Microcoleus", 
                 "Non-Epithemia Diatoms", "Nostoc", "Oscillatoria",
                 "Other Coccoids", "Rare")
colScale <- scale_color_manual(values = mypal)

#Plot

#List of different species groups:
# c("Anabaena", 
#   "Epithemia Diatoms", 
#   "Geitlerinema")

# c("Microcoleus", 
#   "Non-Epithemia Diatoms", 
#   "Green Algae")

mat_p22 <- ggplot(subset(matsims2022, Species %in% c("Microcoleus", 
                                                     "Non-Epithemia Diatoms", 
                                                     "Green Algae")), 
                         aes(x = model_date, y = Abundance)) +
  geom_line(size = 1.5, aes(color = Species)) +
  geom_line(data = subset(mat_params2_TM[mat_params2_TM$year %in% "2022", ], 
                          Species %in% c("Microcoleus", 
                                         "Non-Epithemia Diatoms", 
                                         "Green Algae")), 
            aes(x = model_date, y = mean, colour = Species),
            linewidth = 4, alpha = .35) +
  scale_y_continuous(breaks=c(seq(0,150,5))) +
  #coord_cartesian(ylim = c(0,70)) +
  labs(x = "Date", y = "Relative Abundance (%)", title = "2022 Predictions") +
  colScale


#####################################
#Historical Predictions, 2023
#####################################
#Pull out community abundances and demographics 
abundances <- x[["n"]][,,14] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 13 #number of weeks in 2023
n <- array(NA, dim = c(runs, 11, time))

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
  n[z,,1] <- abundances[z,]
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
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    # #Remove env drivers
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z],
    #                          Sigma = sigma)
    
    
  }
}

pred2023 <- n

matsims2023 <-  as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:11, exp)) %>%
  dplyr::rename(Anabaena = V1, 
                'Epithemia Diatoms' = V2,
                Geitlerinema = V3,
                'Green Algae' = V4, 
                Leptolyngbya = V5,
                Microcoleus = V6,
                'Non-Epithemia Diatoms' = V7,
                Nostoc = V8,
                Oscillatoria = V9,
                'Other Coccoids' = V10,
                Rare = V11) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = c(1:11), names_to = "Species", values_to = "Abundance") %>% #cols are 1:9 if species are ungrouped, 1:3, 5:6 if grouped
  mutate(real_week = time + 26, year = 2023) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))




mat_p23 <- ggplot(subset(matsims2023, Species %in% c("Microcoleus", 
                                                     "Non-Epithemia Diatoms", 
                                                     "Green Algae")), 
                  aes(x = model_date, y = Abundance)) +
  geom_line(size = 1.5, aes(color = Species)) +
  geom_line(data = subset(mat_params2_TM[mat_params2_TM$year %in% "2023", ], 
                          Species %in% c("Microcoleus", 
                                         "Non-Epithemia Diatoms", 
                                         "Green Algae")), 
            aes(x = model_date, y = mean, colour = Species),
            linewidth = 4, alpha = .35) +
  scale_y_continuous(breaks=c(seq(0,150,5))) +
  #coord_cartesian(ylim = c(0,70)) +
  labs(x = "Date", y = "Relative Abundance (%)", title = "2023 Predictions") +
  colScale

#####################################
#Historical Predictions, 2024
#####################################
#Pull out community abundances and demographics 
abundances <- x[["n"]][,,27] #iterations, species #, time

#inputs
runs <- nrow(abundances)
time <- 15 #number of weeks in 2023
n <- array(NA, dim = c(runs, 11, time))

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
  n[z,,1] <- abundances[z,]
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
    n[z,,t] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[z,,t-1] + nTheta*nitrate[t-1]+
                               pTheta*phos[t-1] + aTheta*amon[t-1] + dTheta*dis[t-1] +
                               tTheta*temp[t-1] + cTheta*cond[t-1] + rTheta*rad[t-1],
                             Sigma = sigma)
    # #Remove env drivers
    # n[t,,z] <- MASS::mvrnorm(n = 1, mu = Alpha + Beta%*%n[t-1,,z],
    #                          Sigma = sigma)
    
  }
}

pred2024 <- n

matsims2024 <-  as.data.frame(t(as.data.frame(apply(n, c(2,3), mean)))) %>% 
  dplyr::mutate(across(1:11, exp)) %>%
  dplyr::rename(Anabaena = V1, 
                'Epithemia Diatoms' = V2,
                Geitlerinema = V3,
                'Green Algae' = V4, 
                Leptolyngbya = V5,
                Microcoleus = V6,
                'Non-Epithemia Diatoms' = V7,
                Nostoc = V8,
                Oscillatoria = V9,
                'Other Coccoids' = V10,
                Rare = V11) %>% 
  mutate(time = 1:time) %>% 
  pivot_longer(cols = c(1:11), names_to = "Species", values_to = "Abundance") %>% #cols are 1:9 if species are ungrouped, 1:3, 5:6 if grouped
  mutate(real_week = time + 26, year = 2024) %>% 
  mutate(model_date = ceiling_date(ymd(paste(year, "01", "01", sep = "-")) + 
                                     (real_week - 1) * 7 - 1, "week", week_start = 7))




mat_p24 <- ggplot(subset(matsims2024, Species %in% c("Microcoleus", 
                                                     "Non-Epithemia Diatoms", 
                                                     "Green Algae")), 
                  aes(x = model_date, y = Abundance)) +
  geom_line(size = 1.5, aes(color = Species)) +
  geom_line(data = subset(mat_params2_TM[mat_params2_TM$year %in% "2024", ], 
                          Species %in% c("Microcoleus", 
                                         "Non-Epithemia Diatoms", 
                                         "Green Algae")), 
            aes(x = model_date, y = mean, colour = Species),
            linewidth = 4, alpha = .35) +
  scale_y_continuous(breaks=c(seq(0,150,5))) +
  #coord_cartesian(ylim = c(0,70)) +
  labs(x = "Date", y = "Relative Abundance (%)", title = "2024 Predictions") +
  colScale


ggarrange(
  mat_p22, mat_p23, mat_p24, labels = c("A", "B", "C"), ncol = 3,
  common.legend = TRUE, legend = "bottom"
)

#Compile model check dataframes into a single full timeseries matrix
predictives_mats <- abind(pred2022, pred2023, pred2024, along = 3)



