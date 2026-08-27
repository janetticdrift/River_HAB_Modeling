#CREATE MODEL FOR RIVER-WIDE DATA

#Gather predicted states of percent cover abundances
algaemodelAll <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_AllVar.rds"))
algaemodelBiotic <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Biotic.rds"))
algaemodelAbiotic <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Abiotic.rds"))
algaemodelAbioticnonut <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_AbioticNoNut.rds"))
algaemodelTrueAbiotic <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_TrueAbiotic.rds"))

#Create a list of all the models
abundances_list <- list(
  All = algaemodelAll,
  Biotic = algaemodelBiotic,
  Abiotic = algaemodelAbiotic,
  Abioticnonut = algaemodelAbioticnonut,
  TrueAbiotic = algaemodelTrueAbiotic
)

#Create an empty list that will get filled by the four unique design matrices
X2TM_list <- list()

for (name in names(abundances_list)) {
  
  model <- abundances_list[[name]]
  
  River_predict1 <- as.data.frame(model) %>% 
    t 
  
  River_predict <- as.data.frame(River_predict1) %>% 
    rownames_to_column(var="ID") %>% 
    tidyr::separate_wider_delim(ID, ".", names = c("name", "time")) %>% 
    group_by(name, time) %>% 
    dplyr::mutate(time = as.numeric(time)) %>% 
    dplyr::summarise(median = median(log(c_across(starts_with("V"))), na.rm = TRUE)) %>% 
    dplyr::mutate(Species = case_when(name == 1 ~ 'Green Algae',
                                      name == 2 ~ 'Microcoleus',
                                      name == 3 ~ 'Anabaena',
                                      name == 4 ~ 'Other N Fixers')) %>% 
    ungroup() %>% 
    dplyr::select(!name) %>% 
    pivot_wider(names_from = Species, values_from = median) %>% 
    arrange(time)
  
  
  #Save cleaned percent cover prediction dataframe for making toxin prediction simulation
  saveRDS(River_predict,
          here::here("data/Outputs for Sims and Model Fits/Predicted States/", paste0("River_predict_", name, ".rds")))
  
  #Save cleaned percent cover latent dataframe inside each unique design matrix
  X2TM_list[[name]] <- cbind(
    intercept = 1,
    River_predict[-c(14:15, 29:30), -1],  #Abundances are NOT log-transformed
    nitrate = stand_nut$nitrate_mg_N_L[-c(14:15, 29:30)],
    phos = stand_nut$oPhos_ug_P_L[-c(14:15, 29:30)],
    ammonium = stand_nut$ammonium_mg_N_L[-c(14:15, 29:30)],
    discharge = discharge$stand_discharge[-c(14:15, 29:30)],
    temp = stand_nut$temp_C[-c(14:15, 29:30)],
    cond = stand_nut$cond_uS_cm[-c(14:15, 29:30)],
    rad = swradiation$stand_rad[-c(14:15, 29:30)]
    )
  
  #Create a separate design matrix only for the TAC concentrations using "All" variables
  if(name == "All"){
    X2TAC <- cbind(
      intercept = 1,
      River_predict[-c(1:2, 14:16, 29:33, 41:45), -1],  #Abundances are NOT log-transformed
      nitrate = stand_nut$nitrate_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)],
      phos = stand_nut$oPhos_ug_P_L[-c(1:2, 14:16, 29:33, 41:45)],
      ammonium = stand_nut$ammonium_mg_N_L[-c(1:2, 14:16, 29:33, 41:45)],
      discharge = discharge$stand_discharge[-c(1:2, 14:16, 29:33, 41:45)],
      temp = stand_nut$temp_C[-c(1:2, 14:16, 29:33, 41:45)],
      cond = stand_nut$cond_uS_cm[-c(1:2, 14:16, 29:33, 41:45)],
      rad = swradiation$stand_rad[-c(1:2, 14:16, 29:33, 41:45)]
    )
  }
  
}

#Combine with other model variables into model lists to use in stan() function
#Use toxins samples from Microcoleus mats
#All Variables
model.atx.riverTM_All <- list("uniqueID" = nrow(anatoxin_data_TM),
                          "is_obs" = anatoxin_data_TM$is_obs, #was this an observed day?
                          "firstdays" = anatoxin_data_TM$firstday,
                          "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson needs as.integer
                          "Nspecies" = as.integer(ncol(River_predict)-1),
                          "X2" = X2TM_list$All,
                          "Npredictors" = ncol(X2TM_list$All)
)

model.atx.riverTM_Biotic <- list("uniqueID" = nrow(anatoxin_data_TM),
                              "is_obs" = anatoxin_data_TM$is_obs, #poisson edit, was this an observed day?
                              "firstdays" = anatoxin_data_TM$firstday,
                              "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                              "Nspecies" = as.integer(ncol(River_predict)-1),
                              "X2" = X2TM_list$Biotic,
                              "Npredictors" = ncol(X2TM_list$Biotic)
)

model.atx.riverTM_Abiotic <- list("uniqueID" = nrow(anatoxin_data_TM),
                                 "is_obs" = anatoxin_data_TM$is_obs, #poisson edit, was this an observed day?
                                 "firstdays" = anatoxin_data_TM$firstday,
                                 "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                                 "Nspecies" = as.integer(ncol(River_predict)-1),
                                 "X2" = X2TM_list$Abiotic,
                                 "Npredictors" = ncol(X2TM_list$Abiotic)
)

model.atx.riverTM_AbioticNoNut <- list("uniqueID" = nrow(anatoxin_data_TM),
                                  "is_obs" = anatoxin_data_TM$is_obs, #poisson edit, was this an observed day?
                                  "firstdays" = anatoxin_data_TM$firstday,
                                  "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                                  "Nspecies" = as.integer(ncol(River_predict)-1),
                                  "X2" = X2TM_list$Abioticnonut,
                                  "Npredictors" = ncol(X2TM_list$Abioticnonut)
)

model.atx.riverTM_TrueAbiotic <- list("uniqueID" = nrow(anatoxin_data_TM),
                                       "is_obs" = anatoxin_data_TM$is_obs, #poisson edit, was this an observed day?
                                       "firstdays" = anatoxin_data_TM$firstday,
                                       "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                                       "Nspecies" = as.integer(ncol(River_predict)-1),
                                       "X2" = X2TM_list$TrueAbiotic,
                                       "Npredictors" = ncol(X2TM_list$TrueAbiotic)
)

#Use toxins samples from Anabaena mats
model.atx.riverTAC_All <- list("uniqueID" = nrow(anatoxin_data_TAC),
                           "is_obs" = anatoxin_data_TAC$is_obs, #poisson edit, was this an observed day?
                           "firstdays" = anatoxin_data_TAC$firstday,
                           "Toxins" = as.integer(anatoxin_data_TAC$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                           "Nspecies" = as.integer(ncol(River_predict)-1),
                           "X2" = X2TAC,
                           "Npredictors" = ncol(X2TAC)
)