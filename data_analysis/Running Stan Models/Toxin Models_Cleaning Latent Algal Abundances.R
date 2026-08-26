#CREATE MODEL FOR RIVER-WIDE DATA

#Gather predicted states of percent cover abundances
rivermodelAll <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_AllVar.rds"))
rivermodelBiotic <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Biotic.rds"))
rivermodelAbiotic <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_Abiotic.rds"))
rivermodelAbioticnonut <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_AbioticNoNut.rds"))
rivermodelTrueAbiotic <- readRDS(here::here("data/Outputs for Sims and Model Fits/Predicted States/Riverwide_Pred_TrueAbiotic.rds"))

#Create a list of all the models
abundances_list <- list(
  All = rivermodelAll,
  Biotic = rivermodelBiotic,
  Abiotic = rivermodelAbiotic,
  Abioticnonut = rivermodelAbioticnonut,
  TrueAbiotic = rivermodelTrueAbiotic
)

#Create an empty list that will get filled by the four unique design matrices
X2TM_list <- list()

for (name in names(abundances_list)) {
  
  model <- abundances_list[[name]]
  
  River_latent1 <- as.data.frame(model) %>% 
    dplyr::select(matches("n\\[")) %>% 
    t 
  
  River_latent <- as.data.frame(River_latent1) %>% 
    rownames_to_column(var="ID") %>% 
    tidyr::separate_wider_delim(ID, ".", names = c("chain", "group")) %>% 
    dplyr::select(-chain) %>% 
    group_by(group) %>% 
    dplyr::summarise(median = median(c_across(starts_with("V")), na.rm = TRUE)) %>% 
    dplyr::mutate(Species = case_when(grepl("[1,", group, fixed=TRUE) ~ 'Green Algae',
                                      grepl("[2,", group, fixed=TRUE) ~ 'Microcoleus',
                                      grepl("[3,", group, fixed=TRUE) ~ 'Anabaena',
                                      grepl("[4,", group, fixed=TRUE) ~ 'Other N Fixers')) %>% 
    mutate(time = as.numeric(str_extract_all(group, "[0-9]+", simplify = T)[,2])) %>% 
    dplyr::select(-group) %>% 
    pivot_wider(names_from = Species, values_from = median) %>% 
    arrange(time)
  
  
  #Save cleaned percent cover latent dataframe for making toxin prediction simulation
  saveRDS(River_latent,
          here::here("data/Outputs for Sims and Model Fits/Latent States/", paste0("River_latent_", name, ".rds")))
  
  #Save cleaned percent cover latent dataframe inside each unique design matrix
  X2TM_list[[name]] <- cbind(
    intercept = 1,
    River_latent[-c(14:15, 29:30), -1],  #Abundances are log-transformed
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
      River_latent[-c(1:2, 14:16, 29:33, 41:45), -1],  #Abundances are log-transformed
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
                          "is_obs" = anatoxin_data_TM$is_obs, #poisson edit, was this an observed day?
                          "firstdays" = anatoxin_data_TM$firstday,
                          "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                          "Nspecies" = as.integer(ncol(River_latent)-1),
                          "X2" = X2TM_list$All,
                          "Npredictors" = ncol(X2TM_list$All)
)

model.atx.riverTM_Biotic <- list("uniqueID" = nrow(anatoxin_data_TM),
                              "is_obs" = anatoxin_data_TM$is_obs, #poisson edit, was this an observed day?
                              "firstdays" = anatoxin_data_TM$firstday,
                              "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                              "Nspecies" = as.integer(ncol(River_latent)-1),
                              "X2" = X2TM_list$Biotic,
                              "Npredictors" = ncol(X2TM_list$Biotic)
)

model.atx.riverTM_Abiotic <- list("uniqueID" = nrow(anatoxin_data_TM),
                                 "is_obs" = anatoxin_data_TM$is_obs, #poisson edit, was this an observed day?
                                 "firstdays" = anatoxin_data_TM$firstday,
                                 "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                                 "Nspecies" = as.integer(ncol(River_latent)-1),
                                 "X2" = X2TM_list$Abiotic,
                                 "Npredictors" = ncol(X2TM_list$Abiotic)
)

model.atx.riverTM_AbioticNoNut <- list("uniqueID" = nrow(anatoxin_data_TM),
                                  "is_obs" = anatoxin_data_TM$is_obs, #poisson edit, was this an observed day?
                                  "firstdays" = anatoxin_data_TM$firstday,
                                  "Toxins" = as.integer(anatoxin_data_TM$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                                  "Nspecies" = as.integer(ncol(River_latent)-1),
                                  "X2" = X2TM_list$Abioticnonut,
                                  "Npredictors" = ncol(X2TM_list$Abioticnonut)
)

#Use toxins samples from Anabaena mats
model.atx.riverTAC_All <- list("uniqueID" = nrow(anatoxin_data_TAC),
                           "is_obs" = anatoxin_data_TAC$is_obs, #poisson edit, was this an observed day?
                           "firstdays" = anatoxin_data_TAC$firstday,
                           "Toxins" = as.integer(anatoxin_data_TAC$ATX_all_ug_afdm_g), #poisson edit needs as.integer
                           "Nspecies" = as.integer(ncol(River_latent)-1),
                           "X2" = X2TAC,
                           "Npredictors" = ncol(X2TAC)
)